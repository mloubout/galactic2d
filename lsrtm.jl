using JUDI, Serialization, SlimOptim, Random, LinearAlgebra, MECurvelets
using SlimPlotting

data_path = "/data/galactic2D/"
ver = "l18_b25_mutewater"
im_path = "images/splsrtm/$(ver)"
res_path = "$(data_path)/splsrtm/$(ver)"
mkpath(im_path)
mkpath(res_path)

linenum = 18

include("utils.jl")

# Model and origin
model, origin  = read_model("$(data_path)W22GAL_LINE$(linenum)_FTPreSTM_MigVel.segy"; d=12.5,
                             vals="$(data_path)galactic-CP00054-Vp-0524.sgy", density=true);

wb = find_water_bottom(model.m, 1.6^(-2))
Tm = judiTopmute(model.n, wb, 0)

# OOC data
data, q = load_slice(linenum, "$(data_path)GAL_FullArray_FFSig_Ver2_68pt5ms_0pt5ms.sgy")

C = joMECurvelet2D(model.n...; DDT=Float32, RDT=Float32, zero_finest=true)

# Setup linearized bregman, one datapath
niter = 50
batchsize = div(data.nsrc, 25)
g_scale = 0

# Setup operators
opt = Options(subsampling_factor=8, IC="isic", free_surface=true, limit_m=true, space_order=16)

rand_pos = Set(randperm(q.nsrc))
popn = () -> [pop!(rand_pos) for i=1:batchsize]

function obj(x)
    flush(stdout)
    dm = PhysicalParameter(x, model.n, model.d, model.o)
    if length(rand_pos) > batchsize
        inds = popn()
    else
        union!(rand_pos, Set(randperm(q.nsrc)))
        inds = popn()
    end

    dobs, qf = get_subset(data, q, origin, inds, 5f0, 40f0)
    DG = judiTimeGain(dobs.geometry, 1) # t^2
    Ml = judiDataMute(qf.geometry, dobs.geometry)

    # Imaging
    M = judiModeling(model, qf.geometry, dobs.geometry; options=opt)
    J = judiJacobian(M, qf)
    I = inv(judiIllumination(J; mode="v", k=1))

    dsyn = J*dm
    residual = Ml*(dsyn - DG*dobs)
    # grad
    G = Tm*I*(J'*Ml'*residual)
    g_scale == 0 && (global g_scale = .05f0/norm(G, Inf))
    G .*= g_scale
    gfull = PhysicalParameter(model.n, model.d, model.o, zeros(Float32, model.n))
    gfull .+= G
    return .5f0*norm(residual)^2, gfull[:]
end

cmap = seiscm(:seismic)
perc = 99

function callback(sol)
    serialize("$(res_path)/splsrtm_line18.bin", (x=sol.x, phi=sol.ϕ_trace, z=sol.z))
    i = length(sol.ϕ_trace)-1
    xplot = reshape(sol.x, model.n)
    zplot = reshape(C'*sol.z, model.n)
    figure(figsize=(30, 12));
    plot_simage(xplot', model.d; cmap=cmap, d_scale=0, perc=perc,  new_fig=false, name="iter $(i)", save="$(im_path)/lsrtm_x_$(i)")
    figure(figsize=(30, 12));
    plot_simage(zplot', model.d; cmap=cmap, d_scale=0, perc=perc,  new_fig=false, name="iter $(i)", save="$(im_path)/lsrtm_z_$(i)")
    close("all")
end


dm0 = zeros(Float32, prod(model.n))

# Bregman
bregopt = bregman_options(maxIter=niter, verbose=2, quantile=.95, alpha=1, TD=C,
                          antichatter=false, spg=true, reset_lambda=10)
solb = bregman(obj, dm0, bregopt; callback=callback);
