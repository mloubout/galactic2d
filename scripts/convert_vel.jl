using SegyIO, SlimPlotting, Interpolations

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
        depth[t, :] .= sum(velocity[1:t, :], dims=1)[1, :] .* dt ./ 4
    end
    # depth .*= .5f0

    # Project on regular grid
    maxD = (div(maximum(depth), 12.5) + 1) * 12.5
    @show maxD
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


vp  = segy_read("$(data_path)W22GAL_LINE18_FTPreSTM_MigVel.segy");

# Let's check headers for grid information
X = Float32.(get_header(vp, "CDPX"))
Y = Float32.(get_header(vp, "CDPY"))
dx = (diff(Y)[1]^2 + diff(X)[1]^2)^(.5)

origin = (minimum(X), minimum(Y))

vp_depth, dz = project_velocity(vp)


plot_velocity(vp_depth', (dx, dz); cbar=true, vmax=4500, cmap="gist_ncar")