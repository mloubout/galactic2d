using SegyIO, Serialization, SlimPlotting, JUDI, Interpolations, ImageFiltering

close("all")

include("utils.jl")

linenum = 18
data_path = "/data/galactic2D/"
figpath = "./images/rtm-2405"

model, origin  = read_model("$(data_path)W22GAL_LINE$(linenum)_FTPreSTM_MigVel.segy")
segyvel = "$(data_path)galactic-CP00054-Vp-0524.sgy"

rtm_start = "$(data_path)rtm_line_18_0524_xwi_start_125.bin"
rtm_xwi_125 = "$(data_path)rtm_line_18_0524_xwi_54_125.bin"

cmap = seiscm(:seismic)
vcmap = "cet_rainbow4"
perc = 98
vp_max = 5.0

if isfile("datanorms")
    normsd = deserialize("datanorms")
    norms = normsd.norms
else
    data, q = load_slice(linenum, "/data/galactic2D/GAL_FullArray_FFSig_Ver2_68pt5ms_0pt5ms.sgy")

    norms = [norm(get_data(data[i])) for i=1:data.nsrc]
    itp = LinearInterpolation(range(0, 1, length=data.nsrc), norms, extrapolation_bc=Line())
    norms = itp(range(0, 1, length=model.n[1]))
    serialize("datanorms", (norms=norms,))
end

norms = imfilter(norms, Kernel.gaussian((5,)))

names = ("start", "xwi-54-125")
rtms = (rtm_start, rtm_xwi_125)

for (resname, name) in zip(rtms, names)
    !isfile(resname) && continue
    res = deserialize(resname)
    rtm, Ilu, Ilv, spacing = res.rtm, res.Ilu, res.Ilv, reverse(res.spacing)
    Il = Ilu .* Ilv

    vp = res.m.^(-.5)
    wb = find_water_bottom(res.m, 1.6^(-2))
    Tm = judiTopmute(model.n, wb, 0)

    rtm = Tm * rtm

    nz = size(vp, 2)
    lasti = size(vp, 1)
    startx = 1

    # Illum
    @show extrema(Il[startx:lasti+startx-1, 1:nz])
    figure(figsize=(30, 12));
    plot_velocity(log.(1 .+ Il[startx:lasti+startx-1, 1:nz])', spacing; cmap="gist_ncar", d_scale=0, name="", cbar=true, new_fig=false)
    savefig("$(figpath)/illum_$(name)", bbox_inches = "tight", dpi = 150)

    # RTMs
    figure(figsize=(30, 12));
    plot_simage((rtm[startx:lasti+startx-1, 1:nz])', spacing;
                 cmap=cmap, d_scale=0, perc=perc,  new_fig=false, name="")
    savefig("$(figpath)/RTM_$(name)_novp_noillum", bbox_inches = "tight", dpi = 150)

    figure(figsize=(30, 12));
    plot_simage((rtm[startx:lasti+startx-1, 1:nz])' ./ (Il[startx:lasti+startx-1, 1:nz]' .+ 1e-6), spacing;
                 cmap=cmap, d_scale=0, perc=perc,  new_fig=false, name="")
    savefig("$(figpath)/RTM_$(name)_novp", bbox_inches = "tight", dpi = 150)
    plot_velocity(vp[1:lasti, 1:nz]', spacing; cmap=vcmap, perc=perc, new_fig=false, alpha=.5, vmax=vp_max, name="")
    savefig("$(figpath)/RTM_$(name)", bbox_inches = "tight", dpi = 150)
    close("all")

    dvp = 0 * vp
    dvp[1:lasti, 2:nz] .= diff(vp[1:lasti, 1:nz], dims=2) ./ vp[1:lasti, 2:nz]
    figure(figsize=(30, 12));
    plot_simage(dvp[1:lasti, 1:nz]', spacing;
                cmap="cet_CET_L1", d_scale=0, perc=perc,  new_fig=false, name="", alpha=.5)
    savefig("$(figpath)/dvp$(name)_novp", bbox_inches = "tight", dpi = 150)
    figure(figsize=(30, 12));
    plot_simage((rtm[startx:lasti+startx-1, 1:nz])' ./ (Il[startx:lasti+startx-1, 1:nz]' .+ 1e-6), spacing;
                cmap=seiscm(:seismic), d_scale=0, perc=99, new_fig=false, name="")
    plot_simage(dvp[1:lasti, 1:nz]', spacing;
                cmap="cet_CET_L1", d_scale=0, perc=perc,  new_fig=false, name="", alpha=.5)
    savefig("$(figpath)/dvp$(name)_overlap", bbox_inches = "tight", dpi = 150)

    # block = deepcopy(segy_read(segyvel))
    # nz, nx = size(block.data)
    # block.data[:, :] = rtm[1:nx, 1:nz]'
    # segy_write("$(figpath)/rtm_$(name)vel_2405.segy", block)

end


# SSO

plot_sso = false

if plot_sso

sso = deserialize("$(data_path)sso_line_18_0524.bin")
lasti_sso = size(sso.sso, 2)
nh = size(sso.sso, 1) ÷ 2

figure(figsize=(30, 12));
plot_simage(sso.sso[nh, 1:lasti_sso, :]' ./ (rtm_cshot.Il[1:lasti_sso, 1:nz] .+ 1e-6)', sso.spacing[2:end];
            cmap=cmap, d_scale=0, perc=perc, new_fig=false)
savefig("$(figpath)/SSO_chsot_novp", bbox_inches = "tight", dpi = 150)
plot_velocity(vp_cshot[1:nz, 1:lasti_sso], rtm_start.spacing; cmap=vcmap, perc=perc, new_fig=false, alpha=.5, vmax=vp_max)
savefig("$(figpath)/SSO_chsot", bbox_inches = "tight", dpi = 150)


####################################################
####################################################
####################################################
# Offset gather

ix = Int(div(70000, 12.5)) #div(lasti_sso, 3)
iz = Int(div(4000, 12.5)) #2*div(nz, 3)
h0 = sso.offsets[1]
d = model.d

rtm = sso.sso[nh, 1:lasti_sso, :] ./ (rtm_cshot.Il[1:lasti_sso, 1:nz] .+ 1e-6)
n = size(rtm)

ogz = sso.sso[:, :, iz] ./ (rtm_cshot.Il[1:lasti_sso, iz]' .+ 1e-6)
ogx = sso.sso[:, ix, :] ./ (rtm_cshot.Il[ix, 1:nz]' .+ 1e-6)

fig, axs = subplots(2, 2, figsize=(30, 12), gridspec_kw=Dict(:height_ratios=>[1, 3], :width_ratios=>[3, 1]))
sca(axs[1, 1])
plot_simage(ogz, d; new_fig=false, name="", o=(h0, 0), labels=(:X, :Offset), cmap=seiscm(:seismic), d_scale=0, perc=perc)
hlines(y=0, colors=:b, xmin=0, xmax=(n[1]-1)*d[1], linewidth=1)
vlines(x=(ix-1)*d[1], colors=:b, ymin=h0, ymax=-h0, linewidth=1)
axs[1, 1].set_xticks([])
axs[2, 2].set_xlabel("")
sca(axs[2, 1])
plot_simage(rtm', d; new_fig=false, name="", cmap=seiscm(:seismic), d_scale=0, perc=perc)
vlines(x=(ix-1)*d[1], colors=:b, ymin=0, ymax=(n[2]-1)*d[2], linewidth=1)
hlines(y=(iz-1)*d[2], colors=:b, xmin=0, xmax=(n[1]-1)*d[1], linewidth=1)
sca(axs[2, 2])
plot_simage(ogx', d; new_fig=false, name="", o=(0, h0), labels=(:Offset, :Depth), cmap=seiscm(:seismic), d_scale=0, perc=perc)
vlines(x=0, colors=:b, ymin=0, ymax=(n[2]-1)*d[2], linewidth=1)
hlines(y=(iz-1)*d[2], colors=:b, xmin=h0, xmax=-h0, linewidth=1)
axs[2, 2].set_yticks([])
axs[2, 2].set_ylabel("")
axs[1, 2].set_visible(false)
subplots_adjust(wspace=0, hspace=0)
savefig("$(figpath)/SSO_chsot_ofs", dpi=150, bbox_inches=:tight)

end

