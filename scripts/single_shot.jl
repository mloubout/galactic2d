using SegyIO, JUDI, SlimPlotting, Interpolations, Serialization, SourceEstimation, LinearAlgebra, JLD2
close("all")
data_path = "/data/galactic2D/"
linenum = 18

include("utils.jl")

model, origin  = read_model("$(data_path)W22GAL_LINE$(linenum)_FTPreSTM_MigVel.segy"; d=12.5,
                            vals="$(data_path)galactic-CP00054-Vp-0524.sgy")

figure(figsize=(20, 10))
plot_velocity((model.m').^(-.5); new_fig=false, cbar=true, vmax=4.500, cmap="gist_ncar", save=true)

data, q = load_slice(linenum, "/data/galactic2D/GAL_FullArray_FFSig_Ver2_68pt5ms_0pt5ms.sgy"; t=7500f0)


function synth_data(idx, q, d_obs, model)
    # Get current shot and project geometry
    newshot, newq = get_subset(data, q, origin, idx, 3f0, 40f0)
    # Setup operators with projected objects
    opt = Options(limit_m=true, free_surface=true, space_order=16)  # ~40 GB of memory per source without subsampling
    M = judiModeling(model, newq.geometry, newshot.geometry; options=opt)
    # Inversion
    d_syn = M*newq
    return d_syn
end

idx = [2500]
# newshot, newq = get_subset(data, q, origin, idx, 3f0, 40f0)

# d_syn = synth_data(idx, q, data, model)

# figure(figsize=(10, 10))
# subplot(121)
# plot_sdata(get_data(data[400]); new_fig=false, cmap="seismic", cbar=true, name="400")
# subplot(122)
# plot_sdata(get_data(data[2500]); new_fig=false, cmap="seismic", cbar=true, name="2500")
# savefig("data2_debug", bbox_inches = "tight", dpi = 150)

function image_shot(idx, d_obs, q, model, origin)
    # Get current shot and project geometry
    newshot, newq = get_subset(data, q, origin, idx, 3f0, 40f0)
    DG = judiTimeGain(newshot.geometry, 1) # t^2
    Ml = judiDataMute(newq.geometry, newshot.geometry)
    newshot = Ml*DG*newshot
    # Normalize due to energy imbalances
    newshot ./= norm(newshot, 2)

    # Setup operators with projected objects
    opt = Options(subsampling_factor=6, IC="isic", limit_m=true,
                  free_surface=false, space_order=16)  # ~40 GB of memory per source without subsampling
    M = judiModeling(model, newq.geometry, newshot.geometry; options=opt)

    # Imaging
    J = judiJacobian(M, newq)
    I = judiIllumination(J; mode="uv")
    rtm = J'*newshot
    
    return rtm, I.illums["v"], I.illums["u"]
end


rtm, Iv, Iu = image_shot(idx, data, q, model, origin)


figure(figsize=(30, 12));
plot_velocity(Iu'; cmap="jet", d_scale=0, name="", cbar=true, new_fig=false, vmax=10)
savefig("Ilu_debug", bbox_inches = "tight", dpi = 150)


figure(figsize=(30, 12));
plot_velocity(Iv'; cmap="jet", d_scale=0, name="", cbar=true, new_fig=false, vmax=1e3)
savefig("Ilv_debug", bbox_inches = "tight", dpi = 150)


figure(figsize=(30, 12));
plot_velocity(model.m'.^(-.5); cmap="jet", d_scale=0, name="", cbar=true, new_fig=false, vmax=4.5)
savefig("vp", bbox_inches = "tight", dpi = 150)


# idx = 100

# newshot = get_data(data[idx]; rel_origin=origin, project="2d")
# F = judiFilter(newshot, fmin, fmax)
# DG = judiTimeGain(newshot.geometry, 1) # t^2
# newshot = DG * (F * newshot)

# newq = get_data(q[idx]; rel_origin=origin, project="2d")
# F = judiFilter(newq, fmin, fmax)
# newq = F * newq

# fwi_op = Options(limit_m=true, IC="fwi", free_surface=true, space_order=8, subsampling_factor=6, buffer_size=1000f0)
# I = inv(judiIllumination(model; mode="uv", k=2))
# g_const = 0


# F0 = judiModeling(model, newq.geometry, newshot.geometry; options=fwi_op)
# dsyn = F0 * newq

# dsynn = dsyn / norm(dsyn, Inf)
# newshotn = newshot / norm(newshot, Inf)

# dphi = dsynn - newshotn

# figure(figsize=(10, 10))
# subplot(131)
# plot_sdata(newshotn; new_fig=false, cmap="seismic", cbar=true, name="Obs")
# subplot(132)
# plot_sdata(dsynn; new_fig=false, cmap="seismic", cbar=true, name="Syn")
# subplot(133)
# plot_sdata(dphi; new_fig=false, cmap="seismic",cbar=true, save="images/residual", name="Residual")

# figure(figsize=(10, 10))
# subplot(131)
# plot_sdata(newshotn.data[1][1:500, 1:100], (2f0, 12.5f0); new_fig=false, cmap="seismic", cbar=true, name="Obs")
# subplot(132)
# plot_sdata(dsynn.data[1][1:500, 1:100], (2f0, 12.5f0); new_fig=false, cmap="seismic", cbar=true, name="Syn")
# subplot(133)
# plot_sdata(dphi.data[1][1:500, 1:100], (2f0, 12.5f0); new_fig=false, cmap="seismic",cbar=true, save="images/residual_trunc", name="Residual")


# taxis = 0:newshotn.geometry.dt[1]:newshotn.geometry.t[1]
# figure(figsize=(20, 10))
# subplot(131)
# plot(taxis[1:200], newshotn.data[1][1:200, 1], "k", label="Obs")
# plot(taxis[1:200], dsynn.data[1][1:200, 1], "r", label="Syn")
# legend()
# subplot(132)
# plot(taxis[1:200], newshotn.data[1][1:200, 100], "k", label="Obs")
# plot(taxis[1:200], dsynn.data[1][1:200, 100], "r", label="Syn")
# legend()
# subplot(133)
# plot(taxis[1:200], newshotn.data[1][1:200, 200], "k", label="Obs")
# plot(taxis[1:200], dsynn.data[1][1:200, 200], "r", label="Syn")
# legend()
# savefig("images/residual_trace.png")


# figure()
# plot(newq.data[1][1:200], label="qfilt")
# plot(q.data[1][1:200], label="q")
# legend()
# savefig("images/wavelet")