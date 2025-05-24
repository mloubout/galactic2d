using SegyIO, JUDI, SlimPlotting, PythonCall, PythonPlot
using Interpolations, Serialization, JLD2, ArgParse
using LinearAlgebra, Random, ArgParse


s = ArgParseSettings()
@add_arg_table s begin
    "--case"
        help = "experiement number"
        arg_type = Int
        default = 80
    "--gardner"
        help = "an option without argument, i.e. a flag"
        action = :store_true
end
args = parse_args(ARGS, s)

casenum = "0$(args["case"])"
gardner = args["gardner"]


fmax = 40

plotclose("all")
data_path = "/data/galactic2D/"
figpath = "./images/vti-$(casenum)"
!isdir(figpath) && mkdir(figpath)

linenum = 18

include("utils.jl")

data, q = load_slice(linenum, "/data/galactic2D/GAL_FullArray_FFSig_Ver2_68pt5ms_0pt5ms.sgy"; t=8000f0)

function image_shot(idx, d_obs, q, model, origin; fmin=5f0, fmax=40f0)
    # Get current shot and project geometry
    newshot, newq = get_subset(data, q, origin, idx, fmin, fmax)
    DG = judiTimeGain(newshot.geometry, 1) # t^2
    Ml = judiDataMute(newq.geometry, newshot.geometry)
    newshot = Ml*DG*newshot

    # Setup operators with projected objects
    opt = Options(subsampling_factor=6, IC="isic", limit_m=true,
                  free_surface=true, space_order=16)  # ~40 GB of memory per source without subsampling
    M = judiModeling(model, newq.geometry, newshot.geometry; options=opt)

    # Imaging
    J = judiJacobian(M, newq)
    I = judiIllumination(J; mode="uv")
    rtm = J'*newshot
    
    return rtm, I
end


f32_read(x) = convert(Matrix{Float32}, transpose(segy_read(x).data))

vp = f32_read("$(data_path)/galactic-$(casenum)-Models/galactic-$(casenum)-FINAL-Vp.sgy") ./ 1f3
eps = f32_read("$(data_path)/galactic-$(casenum)-Models/galactic-$(casenum)-TrueEpsilon.sgy")
del = f32_read("$(data_path)/galactic-$(casenum)-Models/galactic-$(casenum)-TrueDelta.sgy")
del[del .> eps] .= .99f0 .* eps[del .> eps]

if gardner
    den = f32_read("$(data_path)/galactic-$(casenum)-Models/galactic-$(casenum)-FINAL-Density-with-DenChange.sgy") ./ 1f3
    densname = "InvDensity"
else
    vtr = pyimport("vtrtool")
    den = pyconvert(Matrix{Float32}, vtr.VTRModel("$(data_path)/galactic-$(casenum)-Models/galactic-$(casenum)-FINAL-Density.vtr").arrays[0]) ./ 1f3
    densname = "Gardner"

end

n = size(vp)
d = (12.5f0, 12.5f0)


figure(figsize=(30, 12));
plot_velocity(vp', d; d_scale=0, new_fig=false, name="VP $(casenum)", save="$(figpath)/vp")
plotclose("all")

# Actual origin from original data
originm, o = get_orig(vp, "$(data_path)W22GAL_LINE$(linenum)_FTPreSTM_MigVel.segy")
@show o, originm

model = Model(n, d, (o, 0f0), vp.^(-2f0); rho=den, epsilon=eps, delta=del)

# model = resample(model, 2)
# @show model.d, model.n
wb = find_water_bottom(model.m, 1.6^(-2))
Tm = judiTopmute(model.n, wb, 0)

rtm = similar(model.m)
fill!(rtm, 0)
Ilu = similar(model.m)
fill!(Ilu, 0.)
Ilv = similar(model.m)
fill!(Ilv, 0.)

filename = "$(data_path)/galactic-$(casenum)-Models/rtm_$(densname)_125_$(fmax).bin"
batch_size = 20

if isfile(filename)
    rtm, _, dloc, oloc, Ilu, Ilv, iter = deserialize(filename)
    rtm = PhysicalParameter(rtm, dloc, oloc)
    Ilu = PhysicalParameter(Ilu, dloc, oloc)
    Ilv = PhysicalParameter(Ilv, dloc, oloc)
    iter += batch_size
else
    iter = 1
end

# idxsrc = vcat(collect.(i:100:data.nsrc for i in 1:100)...)
Random.seed!(data.nsrc)
idxsrc = randperm(data.nsrc)

for it=iter:batch_size:data.nsrc
    i = idxsrc[it:(it+batch_size-1)]
    flush(stdout)
    t1 = @elapsed begin
        rtm_loc, Iloc = image_shot(i, data, q, model, originm; fmin=5f0, fmax=Float32(fmax))
        rtm_loc = Tm * rtm_loc
        if any(isnan.(rtm_loc.data))
            @info "$i has NaN"
            continue
        end
        rtm .+= rtm_loc
        Ilu .+= Iloc.illums["u"]
        Ilv .+= Iloc.illums["v"]
        serialize(filename, (rtm=rtm.data, m=model.m.data, spacing=rtm.d, origin=rtm.o,
                                Ilu=Ilu.data, Ilv=Ilv.data, iter=it))
    end
    println("Shot $(it) of $(data.nsrc) in: $(trunc(t1; digits=3)) s")

    Il = Ilu.^(0.5f0) .* Ilv
    rtm_plot = rtm.data ./ (Il.data .+ 1f-6)
    figure(figsize=(30, 12));
    plot_simage(rtm_plot', rtm.d; cmap="cet_CET_L1", d_scale=0, perc=98,  new_fig=false, name="RTM $(casenum) $(densname)", save="$(figpath)/rtm")
    plotclose("all")

end

