using JUDI, Serialization, SlimOptim, Random, Distributed, LinearAlgebra, TimeProbeSeismic
using PyPlot, SlimPlotting
@everywhere using LinearAlgebra

data_path = "/data/galactic2D/"
linenum = 18

include("utils.jl")

# Model and origin
# model0, origin  = read_model("$(data_path)W22GAL_LINE$(linenum)_FTPreSTM_MigVel.segy"; d=12.5,
#                              vals="$(data_path)galactic-CheckShotModel-A063Start.sgy");

# Model and origin
model0, origin  = read_model("$(data_path)W22GAL_LINE$(linenum)_FTPreSTM_MigVel.segy"; d=12.5,
                              vals="$(data_path)galactic-CP00054-Vp-0524.sgy", density=true);

wb = find_water_bottom(model0.m, 1.6^(-2))
# OOC data
data, q = load_slice(linenum, "$(data_path)GAL_FullArray_FFSig_Ver2_68pt5ms_0pt5ms.sgy")

offsets = -1500f0:spacing(model0, 1):1500f0

function image_shot(idx, d_obs, q, model)
    GC.gc(true)
    # Get current shot and project geometry
    newshot, newq = get_subset(data, q, origin, idx, 5f0, 40f0; normalize=true)

    DG = judiTimeGain(newshot.geometry, 1) # t^2
    Ml = judiDataMute(newq.geometry, newshot.geometry)
    newshot = Ml*DG*newshot

    # Setup operators with projected objects
    opt = Options(subsampling_factor=6, IC="isic", limit_m=true,
                  free_surface=true, space_order=16)  # ~40 GB of memory per source without subsampling
    M = judiModeling(model, newq.geometry, newshot.geometry; options=opt)
    # Inversion
    J = judiJacobian(M, newq, 64, newshot; offsets=offsets, mode=:QR)
    # Mute water layer
    sso = J'*newshot
    for i=1:size(model, 2)
        sso[:, :, 1:wb[i]] .=0
    end

    return sso
end


function plot_current(sso)
    vcmap = "cet_rainbow4"
    cmap = "cet_CET_L2"
    perc = 95
    vp_max = 5.0
    vpplot = model0.m.data.^(-.5)

    lasti_sso = size(sso.data, 2)
    nh = size(sso.data, 1) ÷ 2
    nz = size(sso.data, 3) 
    d = spacing(model0)
    n = size(model0)


    ix = div(lasti_sso, 3)
    iz = 2*div(nz, 3)
    h0 = offsets[1]

    figure(figsize=(30, 12));
    plot_simage(sso.data[nh, 1:lasti_sso, :]', d;
                cmap=cmap, d_scale=1.25, perc=perc, new_fig=false)
    savefig("images/sso-norm/SSO_xwi_novp", bbox_inches = "tight", dpi = 150)
    plot_velocity(vpplot[1:lasti_sso, 1:nz]', d; cmap=vcmap, perc=perc,
                  new_fig=false, alpha=.5, vmax=vp_max)
    savefig("images/sso-norm/SSO_xwi", bbox_inches = "tight", dpi = 150)


    rtm = sso.data[nh, 1:lasti_sso, :]
    ogz = sso.data[:, :, iz]
    ogx = sso.data[:, ix, :]

    fig, axs = subplots(2, 2, figsize=(30, 12), gridspec_kw=Dict(:height_ratios=>[1, 3], :width_ratios=>[3, 1]))
    sca(axs[1, 1])
    plot_simage(ogz, d; new_fig=false, name="", o=(h0, 0), labels=(:X, :Offset), cmap=cmap, d_scale=0, perc=perc)
    hlines(y=0, colors=:b, xmin=0, xmax=(n[1]-1)*d[1], linewidth=1)
    vlines(x=(ix-1)*d[1], colors=:b, ymin=h0, ymax=-h0, linewidth=1)
    axs[1, 1].set_xticks([])
    axs[2, 2].set_xlabel("")
    sca(axs[2, 1])
    plot_simage(rtm', d; new_fig=false, name="", cmap=cmap, d_scale=1.25, perc=perc)
    vlines(x=(ix-1)*d[1], colors=:b, ymin=0, ymax=(n[2]-1)*d[2], linewidth=1)
    hlines(y=(iz-1)*d[2], colors=:b, xmin=0, xmax=(n[1]-1)*d[1], linewidth=1)
    sca(axs[2, 2])
    plot_simage(ogx', d; new_fig=false, name="", o=(0, h0), labels=(:Offset, :Depth), cmap=cmap, d_scale=0, perc=perc)
    vlines(x=0, colors=:b, ymin=0, ymax=(n[2]-1)*d[2], linewidth=1)
    hlines(y=(iz-1)*d[2], colors=:b, xmin=h0, xmax=-h0, linewidth=1)
    axs[2, 2].set_yticks([])
    axs[2, 2].set_ylabel("")
    axs[1, 2].set_visible(false)
    subplots_adjust(wspace=0, hspace=0)
    savefig("images/sso-norm/SSO_xwi_ofs", dpi=150, bbox_inches=:tight)
    close("all")
end

batchsize = try nworkers();catch e 1 end

save_file = "$(data_path)sso_p64_line_18_0524_norm.bin"

if isfile(save_file)
    sso, dloc, oloc, offsets, iter = deserialize(save_file)
    sso = PhysicalParameter(size(sso), dloc, oloc, sso)
    iter += batchsize + 1
    plot_current(sso)
else
    sso = image_shot(1:batchsize, data, q, model0)
    serialize(save_file, (sso=sso.data, spacing=sso.d, origin=sso.o, offsets=offsets, iter=1))
    iter = batchsize + 1
    plot_current(sso)
end

nsrc_sso = data.nsrc

for i=iter:batchsize:nsrc_sso
    idx = i:(i+batchsize-1)
    flush(stdout)
    t1 = @elapsed begin
        sso_loc = image_shot(idx, data, q, model0)
        sso .+= sso_loc
        serialize(save_file, (sso=sso.data, spacing=sso.d, origin=sso.o, offsets=offsets, iter=i))
        plot_current(sso)
    end
    println("Shot ($(idx[1])-$(idx[end])) of $(data.nsrc) in: $(trunc(t1; digits=3)) s")
end
