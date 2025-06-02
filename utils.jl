using LinearAlgebra, Interpolations, JUDI, SegyIO, JLD2, ImageFiltering, Images

export read_model, get_subset, nx, load_slice, get_orig


function get_orig(vp, vp_file, dx=12.5, oind=nothing)
    segvp = segy_read(vp_file)
    X = get_header(segvp, "CDPX")
    Y = get_header(segvp, "CDPY")
    oind = isnothing(oind) ? div(size(vp, 1), 10) : oind
    if Y[2] > Y[1]
        orig = (minimum(X[oind:end]), minimum(Y[oind:end]), 0)
    else
        orig = (minimum(X[oind:end]), maximum(Y[oind:end]), 0)
    end
    return orig, -(oind - 1)*dx
end


function read_model(vp_file::String; start=1, d=12.5, oind=nothing, vals=nothing, density=false)
    T = Float32
    vp  = segy_read(vp_file);
    dx, dz = d, d
    # X axis
    X = get_header(vp, "CDPX")
    Y = get_header(vp, "CDPY")
    # Read file
    if isnothing(vals)
        # Remove last row with zeros
        velocity = Float32.(vp.data)[1:end-1, start:end]
        velocity[velocity .< 1530.0] .= 1530.0
        # Get time axis
        dt = 1e-6 * get_header(vp, "dt")[1] # ms
        nt = get_header(vp, "ns")[1] - 1
        # Time to depth
        depth = similar(velocity)
        depth[1, :] .= 0.
        for t = 2:nt
            depth[t, :] .= depth[t-1, :] + velocity[t-1, :] .* dt
        end
        depth .*= .5f0
        # Somehow need to half it
        depth .*= .5f0

        # Project on regular grid
        maxD = (div(maximum(depth), d) + 1) * d

        # Preallocate the projected velocity array
        new_z = T(0):T(dz):maxD
        nz = length(new_z)
        projected_velocity = zeros(T, nz, size(depth, 2))
        
        for j in 1:size(depth, 2)
            # Interpolate velocity onto regular grid (z, x)
            itp = LinearInterpolation(depth[:, j], velocity[:, j], extrapolation_bc=Line())
            projected_velocity[:, j] = itp(new_z)
        end

        xmin, xmax = extrema(X)
        ymin, ymax = extrema(Y)
        L = sqrt((xmax - xmin)^2 + (ymax - ymin)^2)
        XX = sqrt.((X .- X[1]).^2 .+ (Y .- Y[2]).^2)
    
        newx = 0:dx:L
        grided_vel = zeros(Float32, size(projected_velocity, 1), length(newx))
        for i = 1:size(grided_vel, 1)
            itp = LinearInterpolation(XX, projected_velocity[i, :], extrapolation_bc=Line())
            grided_vel[i, :] = itp(newx)
        end
        grided_vel = 1f-3 .* grided_vel
    else
        # Read values from segy
        vp_in = segy_read(vals)
        grided_vel = 1f-3 .* Float32.(vp_in.data)
    end

    # Origin
    oind = isnothing(oind) ? div(size(grided_vel, 1), 10) : oind
    if Y[2] > Y[1]
        orig = (minimum(X[oind:end]), minimum(Y[oind:end]), 0)
    else
        orig = (minimum(X[oind:end]), maximum(Y[oind:end]), 0)
    end

    if density
        rho = Gardner(grided_vel; vwater=1.599)
        @show extrema(rho)
        @show extrema(grided_vel)
        model = Model(size(grided_vel'), (dx, dz), (-(oind - 1)*dx, 0.), (grided_vel').^(-2), collect(rho'); nb=80)
    else
        model = Model(size(grided_vel'), (dx, dz), (-(oind - 1)*dx, 0.), (grided_vel').^(-2); nb=80)
    end
    return model, orig
end


function resample(model, ratio)
    newn = (size(model) .- 1) .* ratio .+ 1
    newd = spacing(model) ./ 2
    newo = origin(model)
    m = imresize(model.m.data, newn)
    if isa(model.rho, JUDI.PhysicalParameter)
        rho = imresize(model.rho.data, newn)
        return Model(newn, newd, newo, m, rho; nb=nbl(model))
    else
        return Model(newn, newd, newo, m; nb=nbl(model))
    end
end

nx(x) = x ./ norm(x, Inf)

function load_slice(linenum, data_path, segy_key, wavelet_path::String; t=nothing, t0rec=nothing, t0src=-14f0, src_depthkey="SourceDepth")
    # Data
    if isfile("$(data_path)$(linenum)$(segy_key).jld2")
        @load "$(data_path)$(linenum)$(segy_key).jld2" shots
    else
        shots = segy_scan(data_path, segy_key, ["GroupX", "GroupY", "dt", "ns", "RecGroupElevation", src_depthkey]);
        @save "$(data_path)$(linenum)$(segy_key).jld2" shots
    end
    data = judiVector(shots; segy_depth_key="RecGroupElevation", t=t, t0=t0rec)

    # Source
    # Set up wavelet
    src_geometry = Geometry(shots; key = "source", segy_depth_key=src_depthkey, t=t, t0=t0src)
    wavelet = segy_read(wavelet_path)
    dtw = get_header(wavelet, "dt")[1]/1000
    nsw = get_header(wavelet, "ns")[1]
    twavelet = 0:dtw:((nsw-1)*dtw)

    dtd = get_dt(src_geometry, 1)
    newt = 0:dtd:twavelet[end]
    @info "Input wavelet with dt=$(dtw), nt=$(nsw), data sampling=$(dtd)"

    # Hack for scube wavelet
    if contains(wavelet_path, "Vendor")
        i0 = Int(nsw / 2) - 100
        src_data = Float32.(wavelet.data)[i0:end]
        twavelet = 0:dtw:((length(src_data)-1)*dtw)
    else
        src_data = Float32.(wavelet.data)[:]
    end
    
    itp = LinearInterpolation(twavelet, src_data, extrapolation_bc=Line())
    itq = itp(newt)
    ntd = get_nt(src_geometry, 1)
    wavelet_q = zeros(Float32, ntd, 1)
    ntq = min(ntd, length(itq))
    wavelet_q[1:ntq] .= itq[1:ntq]

    q = -diff(judiVector(src_geometry, wavelet_q), dims=1)

    return data, q
end

function get_subset(data, q, orig, idx, f0=3f0, f1=30f0; normalize=false)
    newq = get_data(q[idx]; rel_origin=orig, project="2d")
    newshot = get_data(data[idx]; rel_origin=orig, project="2d")
    if isnothing(f0) || isnothing(f1)
        return newshot, newq
    end
    Fq = judiFilter(newq, f0, f1)
    Fd = judiFilter(newshot, f0, f1)
    newq = Fq * newq
    newshot = Fd * newshot
    if normalize
        newq ./= map(norm, newq)
        newshot ./= map(norm, newshot)
    end
    return newshot, newq
end