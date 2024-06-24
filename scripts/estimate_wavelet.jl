using SegyIO, JUDI, SlimPlotting, Interpolations, Serialization, SourceEstimation

data_path = "/data/galactic2D/"
function project_velocity(vp::SeisBlock)
    T = Float32
    # Remove last row with zeros
    velocity = Float32.(vp.data)[1:end-1, :]
    # Get time axis
    dt = 1e-6 * get_header(vp, "dt")[1] # ms
    nt = get_header(vp, "ns")[1] - 1
    time_axis = Float32.(0:dt:((nt-2)*dt));
    # Time to depth
    depth = similar(velocity)
    depth[1, :] .= 0.
    for t = 2:nt
        depth[t, :] .= depth[t-1, :] + .5 * (velocity[t, :] .+ velocity[t-1, :]) .* dt
    end
    depth .*= .5f0

    # Project on regular grid
    maxD = (div(maximum(depth), 12.5) + 1) * 12.5

    # Preallocate the projected velocity array
    new_z = T(0):12.5f0:maxD
    nz = length(new_z)
    projected_velocity = zeros(T, nz, size(depth, 2))
    
    for j in 1:size(depth, 2)
        # Interpolate velocity onto regular grid (z, x)
        itp = LinearInterpolation(depth[:, j], velocity[:, j], extrapolation_bc=Line())
        projected_velocity[:, j] = itp(new_z)
    end
    
    return collect(projected_velocity'), diff(new_z)[1]
end

vp  = segy_read("$(data_path)W22GAL_LINE01_FTPreSTM_MigVel.segy");

# Let's check headers for grid information
X = Float32.(get_header(vp, "CDPX"))
Y = Float32.(get_header(vp, "CDPY"))
dx = (diff(Y)[1]^2 + diff(X)[1]^2)^(.5)

origin = (minimum(X), minimum(Y))

vp_depth, dz = project_velocity(vp)

figure(figsize=(20, 10))
plot_velocity(vp_depth', (dx, dz); new_fig=false, cbar=true, vmax=4500, cmap="gist_ncar")

model = Model(size(vp_depth), (dx, dz), (0., 0.), (1e-3*vp_depth).^(-2))

shots = segy_scan(data_path, "SHOTS", ["GroupX", "GroupY", "dt", "ns", "SourceX", "SourceY"]);

data = judiVector(shots)

shot10 = get_data(data[10]);
F = judiFilter(shot10, 3f0, 50f0)

figure(figsize=(10, 10))
plot_sdata(F*shot10; cmap="seismic", new_fig=false)

# Set up wavelet
src_geometry = Geometry(shots; key = "source", segy_depth_key = "SourceDepth")
wavelet = ricker_wavelet(src_geometry.t[1], src_geometry.dt[1], 0.015)    # 30 Hz wavelet
q = judiVector(src_geometry, wavelet)

function project_geometry(G::GeometryIC{T}) where T
    new_x = sqrt.((G.xloc[1] .- origin[1]).^2 .+ (G.yloc[1] .- origin[2]).^2)
    new_y = 0 .* new_x
    return Geometry([new_x], [new_y], G.zloc, G.dt, G.nt, G.t)
end

idx = 10
# Get current shot and project geometry
newshot = judiVector(project_geometry(Geometry(data[idx].geometry)), Float32.(data.data[idx][1].data))
newq = judiVector(project_geometry(Geometry(q[idx].geometry)), q[idx].data)
# Setup operators with projected objects
opt = Options(limit_m=true)  # ~40 GB of memory per source without subsampling
M = judiModeling(model, newq.geometry, newshot.geometry; options=opt)
# Inversion
F = judiFilter(newshot, 3f0, 50f0)

d_syn = M*newq
P, qest = estimate_source!(d_syn, newshot, newq; invQ_TaperWidth=200, invQ_lambda1=1f-3, invQ_lambda2=1f-3)
