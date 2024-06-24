using JUDI, Serialization, SlimOptim, Random, Distributed, LinearAlgebra, SlimPlotting
@everywhere using LinearAlgebra

data_path = "/data/galactic2D/"
linenum = 18

include("utils.jl")

# Model and origin
model, origin  = read_model("$(data_path)W22GAL_LINE$(linenum)_FTPreSTM_MigVel.segy"; d=25.0);

# OOC data
data, q = load_slice(linenum, "$(data_path)GAL_FullArray_FFSig_Ver2_68pt5ms_0pt5ms.sgy"; t=5500f0)

# Inversion options
fwi_op = Options(limit_m=true, IC="fwi", free_surface=true, space_order=8, subsampling_factor=6)
g_const = 0

# Phase only misfit (alignment) wich has a normalized misfit as adjoint source
@everywhere function phasefit(dsyn, dobs)
    scale = norm(dsyn[1:80, 1:80], Inf) / norm(dobs[1:80, 1:80], Inf)
    dobs *= scale
    phi = .5f0*norm(dsyn - dobs)^2
    dphi = dsyn - dobs
    return phi, dphi
end

function phasefit(dsyn::judiVector, dobs::judiVector)
    fval = 0
    for i=1:dsyn.nsrc
        scale = norm(dsyn.data[i][1:150, :], Inf) / norm(dobs.data[i][1:150, :], Inf)
        fval += .5f0*norm(dsyn[i] - scale*dobs[i])^2
    end
    return fval
end

"""
Objective function
"""
function objective_function(x, G, batchsize, data, q, fmax)
    model.m .= reshape(x, model.n);
    batchsizeloc = isnothing(G) ? div(batchsize, 10) : batchsize

    # Data and preprocess
    idx = randperm(data.nsrc)[1:batchsizeloc]
    dobs, qf = get_subset(data, q, origin, idx, 3f0, fmax)
    DG = judiTimeGain(dobs.geometry, 1) # t^2
    dobs = DG*dobs

    # fwi function value and gradient
    if isnothing(G)
        # FUnction only, e.g for linesearch
        F0 = judiModeling(model, qf.geometry, dobs.geometry; options=fwi_op)
        dsyn = F0*qf
        fval = phasefit(dsyn, dobs)
    else
        I = judiIllumination(model; mode="v", k=2)
        fval, grad = fwi_objective(model, qf, dobs; options=fwi_op, misfit=phasefit)
        G .= inv(I)*grad
        if g_const == 0
            global g_const = .5f0 ./ norm(G, Inf) 
        end
        G .*= g_const

        figure(figsize=(30,12))
        plot_simage(G'; save="./fwi18/Gfwi.png", cmap="seismic", d_scale=0, perc=98, new_fig=false, cbar=true)

        Iplot = inv(I) * (0f0 .* G .+ 1f0)

        figure(figsize=(30,12))
        plot_velocity(Iplot'; save="./fwi18/Illum.png", cmap="hsv", d_scale=0, vmax=1000*minimum(Iplot), perc=98, new_fig=false, cbar=true)
        close("all")
    end

    return fval / batchsizeloc
end

# InvM
invM = deepcopy(model.m)

# Parameters
batchsize = div(data.nsrc, 10)
freqs = [5f0, 7f0, 9f0, 11f0]
iters = [10, 10, 10, 10]

# Bounds constraints
minx, maxx = (5.0)^(-2), (1.48)^(-2)

function proj_bounds(x)
    out = 1 .* x
    out[x .< minx] .= minx
    out[x .> maxx] .= maxx
    return out
end

# Frequency continuation
for (fi, fmax) in enumerate(freqs)

    function callback(sol)
        iter = length(sol.ϕ_trace)-1
        name = "./fwi18/fwi_$(fi)_$(iter).png"
        current = (sol.x).^(-.5f0)

        figure(figsize=(30,12))
        plot_velocity(current'; save=name, cmap="gist_ncar", d_scale=0, vmax=4.5, new_fig=false, cbar=true)
    end

    # Setup functions
    global g_const = 0
    f(x) = objective_function(x, nothing, batchsize, data, q, fmax)
    g!(g, x) = objective_function(x, g, batchsize, data, q, fmax)
    fg!(g, x) = objective_function(x, g, batchsize, data, q, fmax)
    spg_opt = spg_options(memory=4, maxIter=iters[fi], store_trace=true)

    sol = spg(f, g!, fg!, invM, proj_bounds, spg_opt, callback=callback)
    # Save solution
    serialize("$(data_path)fwiL18_$(fi).bin", (sol=sol, fmax=fmax))
    # Update invM
    invM .= sol.x
end