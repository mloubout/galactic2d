using SegyIO, JUDI, SlimPlotting, Interpolations, Serialization, SourceEstimation, LinearAlgebra, JLD2
close("all")
data_path = "/data/galactic2D/"
linenum = 18

include("utils.jl")

model, origin  = read_model("$(data_path)W22GAL_LINE$(linenum)_FTPreSTM_MigVel.segy"; d=12.5);

figure(figsize=(20, 10))
plot_velocity((model.m').^(-.5); new_fig=false, cbar=true, vmax=4.500, cmap="gist_ncar", save=true)

data, q = load_slice(linenum, "/data/galactic2D/GAL_FullArray_FFSig_Ver2_68pt5ms_0pt5ms.sgy"; t=5000f0)

fmin, fmax = 3f0, 30f0

function synth_data(idx, q, d_obs, model)
    # Get current shot and project geometry
    newq = get_data(q[idx]; rel_origin=origin, project="2d")
    newshot = get_data(d_obs[idx]; rel_origin=origin, project="2d")
    Fq = judiFilter(newq, fmin, fmax)
    # Setup operators with projected objects
    opt = Options(limit_m=true, free_surface=true, space_order=16)  # ~40 GB of memory per source without subsampling
    M = judiModeling(model, newq.geometry, newshot.geometry; options=opt)
    # Inversion
    d_syn = M*(Fq*newq)
    return d_syn
end

idx = 100

newshot = get_data(data[idx]; rel_origin=origin, project="2d")
F = judiFilter(newshot, fmin, fmax)
DG = judiTimeGain(newshot.geometry, 1) # t^2
newshot = DG * (F * newshot)

newq = get_data(q[idx]; rel_origin=origin, project="2d")
F = judiFilter(newq, fmin, fmax)
newq = F * newq

fwi_op = Options(limit_m=true, IC="fwi", free_surface=true, space_order=8, subsampling_factor=6, buffer_size=1000f0)
I = inv(judiIllumination(model; mode="uv", k=2))
g_const = 0


F0 = judiModeling(model, newq.geometry, newshot.geometry; options=fwi_op)
dsyn = F0 * newq

dsynn = dsyn / norm(dsyn, Inf)
newshotn = newshot / norm(newshot, Inf)

dphi = dsynn - newshotn

figure(figsize=(10, 10))
subplot(131)
plot_sdata(newshotn; new_fig=false, cmap="seismic", cbar=true, name="Obs")
subplot(132)
plot_sdata(dsynn; new_fig=false, cmap="seismic", cbar=true, name="Syn")
subplot(133)
plot_sdata(dphi; new_fig=false, cmap="seismic",cbar=true, save="images/residual", name="Residual")

figure(figsize=(10, 10))
subplot(131)
plot_sdata(newshotn.data[1][1:500, 1:100], (2f0, 12.5f0); new_fig=false, cmap="seismic", cbar=true, name="Obs")
subplot(132)
plot_sdata(dsynn.data[1][1:500, 1:100], (2f0, 12.5f0); new_fig=false, cmap="seismic", cbar=true, name="Syn")
subplot(133)
plot_sdata(dphi.data[1][1:500, 1:100], (2f0, 12.5f0); new_fig=false, cmap="seismic",cbar=true, save="images/residual_trunc", name="Residual")


taxis = 0:newshotn.geometry.dt[1]:newshotn.geometry.t[1]
figure(figsize=(20, 10))
subplot(131)
plot(taxis[1:200], newshotn.data[1][1:200, 1], "k", label="Obs")
plot(taxis[1:200], dsynn.data[1][1:200, 1], "r", label="Syn")
legend()
subplot(132)
plot(taxis[1:200], newshotn.data[1][1:200, 100], "k", label="Obs")
plot(taxis[1:200], dsynn.data[1][1:200, 100], "r", label="Syn")
legend()
subplot(133)
plot(taxis[1:200], newshotn.data[1][1:200, 200], "k", label="Obs")
plot(taxis[1:200], dsynn.data[1][1:200, 200], "r", label="Syn")
legend()
savefig("images/residual_trace.png")


figure()
plot(newq.data[1][1:200], label="qfilt")
plot(q.data[1][1:200], label="q")
legend()
savefig("images/wavelet")