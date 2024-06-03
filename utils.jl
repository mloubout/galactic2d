using LinearAlgebra, Interpolations, JUDI, SegyIO, JLD2

export read_model, get_subset, nx, load_slice

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
        origin = (minimum(X[oind:end]), minimum(Y[oind:end]), 0)
    else
        origin = (minimum(X[oind:end]), maximum(Y[oind:end]), 0)
    end

    if density
        rho = Gardner(grided_vel; vwater=1.599)
        @show extrema(rho)
        @show extrema(grided_vel)
        model = Model(size(grided_vel'), (dx, dz), (-(oind - 1)*dx, 0.), (grided_vel').^(-2), collect(rho'); nb=80)
    else
        model = Model(size(grided_vel'), (dx, dz), (-(oind - 1)*dx, 0.), (grided_vel').^(-2); nb=80)
    end
    return model, origin
end

nx(x) = x ./ norm(x, Inf)

function load_slice(linenum, wavelet::String; t=nothing)
    # Data
    if isfile("$(linenum).jld2")
        @load "$(linenum).jld2" shots
    else
        shots = segy_scan(data_path, "W22GAL$(linenum)", ["GroupX", "GroupY", "dt", "ns", "RecGroupElevation", "SourceDepth"]);
        @save "$(linenum).jld2" shots
    end
    data = judiVector(shots; segy_depth_key="RecGroupElevation", t=t)
    
    # Source
    # Set up wavelet
    src_geometry = Geometry(shots; key = "source", segy_depth_key = "SourceDepth", t=t, t0=-14f0)
    wavelet = segy_read(wavelet)
    dtw = get_header(wavelet, "dt")[1]/1000
    nsw = get_header(wavelet, "ns")[1]
    twavelet = 0:dtw:((nsw-1)*dtw)

    dtd = get_dt(src_geometry, 1)
    newt = 0:dtd:twavelet[end]
    
    src_data = Float32.(wavelet.data)[:]
    
    itp = LinearInterpolation(twavelet, src_data, extrapolation_bc=Line())
    itq = itp(newt)
    wavelet_q = zeros(Float32, get_nt(src_geometry, 1), 1)
    wavelet_q[1:length(itq)] .= itq
    
    q = -diff(judiVector(src_geometry, wavelet_q), dims=1)

    return data, q
end

function get_subset(data, q, origin, idx, f0=3f0, f1=30f0; normalize=false)
    newq = get_data(q[idx]; rel_origin=origin, project="2d")
    newshot = get_data(data[idx]; rel_origin=origin, project="2d")
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