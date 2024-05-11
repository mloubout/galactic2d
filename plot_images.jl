using SegyIO, Serialization, SlimPlotting, JUDI

close("all")

include("utils.jl")

linenum = 18
data_path = "/data/galactic2D/"

model, origin  = read_model("$(data_path)W22GAL_LINE$(linenum)_FTPreSTM_MigVel.segy")
vp_start = (model.m.data.^(-.5))
vp_check_segy = segy_read("$(data_path)galactic-CheckShotModel-A063Start.sgy")
vp_cshot = Float32.(vp_check_segy.data)' ./ 1f3

rtm_start = "$(data_path)rtm_line_18.bin"
rtm_cshot = "$(data_path)rtm_line_18_sc.bin"
rtm_cshot_2 = "$(data_path)rtm_line_18_xwi_start.bin"
rtm_xwi = "$(data_path)rtm_line_18_xwi_54.bin"
rtm_cshot_2_125 = "$(data_path)rtm_line_18_xwi_start_125.bin"
rtm_xwi_125 = "$(data_path)rtm_line_18_xwi_54_125.bin"

cmap = "cet_CET_L2" #seiscm(:seismic)
vcmap = "cet_rainbow4"
perc = 95
vp_max = 5.0

names = ("start", "cshot", "cshot-xwi", "xwi-64", "cshot-xwi-125", "xwi-64-125")
rtms = (rtm_start, rtm_cshot, rtm_cshot_2, rtm_xwi, rtm_cshot_2_125, rtm_xwi_125)
velocities = ("$(data_path)galactic-CheckShotModel-A063Start.sgy",
              "$(data_path)galactic-CheckShotModel-A063Start.sgy",
              "$(data_path)galactic-A85-StartVp.sgy",
              "$(data_path)galactic-A85-CP00054-Vp.sgy",
              "$(data_path)galactic-A85-StartVp.sgy",
              "$(data_path)galactic-A85-CP00054-Vp.sgy")

for (resname, name, segyvel) in zip(rtms[6:6], names[6:6], velocities[6:6])
    !isfile(resname) && continue
    res = deserialize(resname)
    rtm, Il, spacing = res.rtm, res.Il, reverse(res.spacing)
    if name == "cshot"
        vp = vp_cshot
    elseif name == "start"
        vp = vp_start
    else
        vp = res.m.^(-.5)
        vp = vp
    end
    nz = size(vp, 2)
    lasti = size(vp_start, 1)
    startx = name == "start" ? 201 : 1
    figure(figsize=(30, 12));
    plot_simage((rtm[startx:lasti+startx-1, 1:nz] ./ (Il[startx:lasti+startx-1, 1:nz] .+ 1e-6))', spacing;
                 cmap=cmap, d_scale=0, perc=perc,  new_fig=false, name="")
    savefig("images/RTM_$(name)_novp", bbox_inches = "tight", dpi = 150)
    plot_velocity(vp[1:lasti, 1:nz]', spacing; cmap=vcmap, perc=perc, new_fig=false, alpha=.5, vmax=vp_max, name="RTM $(name)")
    savefig("images/RTM_$(name)", bbox_inches = "tight", dpi = 150)
    close("all")

    dvp = 0 * vp
    dvp[1:lasti, 2:nz] .= diff(vp[1:lasti, 1:nz], dims=2) ./ vp[1:lasti, 2:nz]
    figure(figsize=(30, 12));
    plot_simage(dvp[1:lasti, 1:nz]', spacing;
                 cmap=cmap, d_scale=0, perc=perc,  new_fig=false, name="")
    savefig("images/dvp$(name)_novp", bbox_inches = "tight", dpi = 150)

    # block = deepcopy(segy_read(segyvel))
    # nz, nx = size(block.data)
    # block.data[:, :] = rtm[1:nx, 1:nz]'
    # segy_write("rtm_$(name)vel.segy", block)

end



# SSO

plot_sso = false

if plot_sso

sso = deserialize("$(data_path)sso_line_18.bin")
lasti_sso = size(sso.sso, 2)
nh = size(sso.sso, 1) ÷ 2

figure(figsize=(30, 12));
plot_simage(sso.sso[nh, 1:lasti_sso, :]' ./ (rtm_cshot.Il[1:lasti_sso, 1:nz] .+ 1e-6)', sso.spacing[2:end];
            cmap=cmap, d_scale=0, perc=perc, new_fig=false)
savefig("images/SSO_chsot_novp", bbox_inches = "tight", dpi = 150)
plot_velocity(vp_cshot[1:nz, 1:lasti_sso], rtm_start.spacing; cmap=vcmap, perc=perc, new_fig=false, alpha=.5, vmax=vp_max)
savefig("images/SSO_chsot", bbox_inches = "tight", dpi = 150)


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
savefig("images/SSO_chsot_ofs", dpi=150, bbox_inches=:tight)

end

