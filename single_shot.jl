using SegyIO, JUDI, SlimPlotting, Interpolations, Serialization, SourceEstimation, LinearAlgebra, JLD2
close("all")
data_path = "/data/galactic2D/"
linenum = 18

include("utils.jl")

model, orig  = read_model("$(data_path)W22GAL_LINE$(linenum)_FTPreSTM_MigVel.segy"; d=12.5,
                            vals="$(data_path)galactic-CP00054-Vp-0524.sgy", density=true)

figure(figsize=(20, 10))
plot_velocity((model.m').^(-.5); new_fig=false, cbar=true, vmax=4.500, cmap="gist_ncar", save=true)

data, q = load_slice(linenum, "/data/galactic2D/GAL_FullArray_FFSig_Ver2_68pt5ms_0pt5ms.sgy")


idx = [3447]
d_obs, a_loc = get_subset(data, q, orig, idx, 5f0, fmax)

sx = 1f-3 * a_loc.geometry.xloc
sx = trunc.(Int, [s[1] for s in sx])

figure(figsize=(12, 8))
subplot(131)
plot_sdata(d_obs[1]; new_fig=false, cmap="seismic", cbar=true, name="Source at $(sx[1][1])km")
subplot(132)
plot_sdata(d_obs[2]; new_fig=false, cmap="seismic", cbar=true, name="Source at $(sx[2][1])km")
subplot(133)
plot_sdata(d_obs[3]; new_fig=false, cmap="seismic", cbar=true, name="Source at $(sx[3][1])km")
tight_layout()
savefig("data-arti.png",  bbox_inches = "tight", dpi = 150)


# model = resample(model, 2)
sfactor = 10
fmax = 60f0

function synth_data(idx, q, d_obs, model)
    # Get current shot and project geometry
    newshot, newq = get_subset(data, q, orig, idx, 5f0, fmax)
    DG = judiTimeGain(newshot.geometry, 1) # t^2
    newshot = DG*newshot

    # Setup operators with projected objects
    opt = Options(limit_m=true, free_surface=true, space_order=16)  # ~40 GB of memory per source without subsampling
    M = judiModeling(model, newq.geometry, newshot.geometry; options=opt)
    # Inversion
    d_syn = M*newq
    return d_syn/norm(d_syn, Inf), newshot/norm(newshot, Inf), newq
end

idx = [3447]
# newshot, newq = get_subset(data, q, orig, idx, 3f0, 40f0)

d_syn, d_obs, newq = synth_data(idx, q, data, model)

figure(figsize=(10, 10))
subplot(131)
plot_sdata(d_obs; new_fig=false, cmap="seismic", cbar=true, name="field")
subplot(132)
plot_sdata(d_syn; new_fig=false, cmap="seismic", cbar=true, name="synth")
subplot(133)
plot_sdata(d_syn-d_obs; new_fig=false, cmap="seismic", cbar=true, name="diff")
savefig("data2_debug", bbox_inches = "tight", dpi = 150)

figure()
plot(d_syn.data[1][1:400, 1], label="syn")
plot(d_obs.data[1][1:400, 1], label="obs")
legend()
savefig("firsttrace")

figure()
plot(newq.data[1][1:400])
savefig("source_70hz")

function image_shot(idx, d_obs, q, model, orig)
    # Get current shot and project geometry
    newshot, newq = get_subset(data, q, orig, idx, 5f0, fmax)
    DG = judiTimeGain(newshot.geometry, 1) # t^2
    Ml = judiDataMute(newq.geometry, newshot.geometry)
    newshot = Ml*DG*newshot

    # Setup operators with projected objects
    opt = Options(subsampling_factor=sfactor, IC="isic", limit_m=true,
                  free_surface=true, space_order=16)  # ~40 GB of memory per source without subsampling
    M = judiModeling(model, newq.geometry, newshot.geometry; options=opt)

    # Imaging
    J = judiJacobian(M, newq)
    I = judiIllumination(J; mode="uv")
    rtm = J'*newshot
    
    return rtm, I.illums["v"], I.illums["u"]
end


rtm, Iv, Iu = image_shot(idx, data, q, model, orig)


figure(figsize=(30, 12));
plot_velocity(Iu'; cmap="jet", d_scale=0, name="", cbar=true, new_fig=false)
savefig("Ilu_debug", bbox_inches = "tight", dpi = 150)


figure(figsize=(30, 12));
plot_velocity(Iv'; cmap="jet", d_scale=0, name="", cbar=true, new_fig=false)
savefig("Ilv_debug", bbox_inches = "tight", dpi = 150)


figure(figsize=(30, 12));
plot_velocity(model.m'.^(-.5); cmap="jet", d_scale=0, name="", cbar=true, new_fig=false)
savefig("vp", bbox_inches = "tight", dpi = 150)

@show extrema(rtm.data)
if any(isnan.(rtm.data))
    figure(figsize=(30, 12));
    plot_simage(isnan.(rtm.data)', spacing(model); cmap="cet_CET_L2", d_scale=0, name="", cbar=true, new_fig=false, perc=99)
    savefig("nanrtm", bbox_inches = "tight", dpi = 150)
end

figure(figsize=(30, 12));
plot_simage(rtm'; cmap="cet_CET_L2", d_scale=0, name="", cbar=true, new_fig=false, perc=99)
savefig("singlertm", bbox_inches = "tight", dpi = 150)



# idx = 100

# newshot = get_data(data[idx]; rel_origin=orig, project="2d")
# F = judiFilter(newshot, fmin, fmax)
# DG = judiTimeGain(newshot.geometry, 1) # t^2
# newshot = DG * (F * newshot)

# newq = get_data(q[idx]; rel_origin=orig, project="2d")
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