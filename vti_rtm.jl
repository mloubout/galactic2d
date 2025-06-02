using SegyIO, JUDI, SlimPlotting, PythonCall, PythonPlot
using Interpolations, Serialization, JLD2, ArgParse, ImageFiltering
using LinearAlgebra, Random, ArgParse

import JUDI: IsoElModel

vtr = pyimport("vtrtool")

s = ArgParseSettings()
@add_arg_table s begin
    "--case"
        help = "experiement number"
        arg_type = Int
        default = 80
    "--gardner"
        help = "Use gardner density"
        action = :store_false
    "--no-multiples"
        help = "Data with multiples"
        action = :store_false
    "--vendor"
        help = "Vendor number"
        arg_type = Int
        default = 1
end
args = parse_args(ARGS, s)

casenum = "0$(args["case"])"
gardner = args["gardner"]
fs = args["no-multiples"]
vendor = args["vendor"]


println("Running case $(casenum) with Gardner=$(gardner), fs=$(fs)")

fmax = 40

plotclose("all")
data_path = "/data/galactic2D/"
if !fs
    shot_path = "$(data_path)demultiple/"
    wavelet = "$(data_path)demultiple/Vendor0$(vendor)-Wavelet.sgy"
    t0src = -396f0
    t0rec = 0f0
    multiples = "_vendor$(vendor)"
    segy_key = "SRME_Vendor0$(vendor)"
    if vendor == 2
        src_depthkey = "SourceSurfaceElevation"
    else
        src_depthkey = "SourceDepth"
    end
else
    shot_path = data_path
    wavelet = "/data/galactic2D/GAL_FullArray_FFSig_Ver2_68pt5ms_0pt5ms.sgy"
    t0src = -14f0
    t0rec = nothing
    multiples = ""
    segy_key = "W22GAL$(linenum)P1002"
    src_depthkey = "SourceDepth"
end

figpath = "./images/vti-$(casenum)"
!isdir(figpath) && mkdir(figpath)

linenum = 18

include("utils.jl")

data, q = load_slice(linenum, shot_path, segy_key, wavelet; t=8000f0, t0src=t0src, t0rec=t0rec, src_depthkey=src_depthkey)

ic_mode(m) = "isic"
ic_mode(::IsoElModel) = "as"

post_process(x::Matrix{T}, model) where T = x
post_process(x::Matrix{T}, ::IsoElModel) where T = convert(Matrix{T}, imfilter(x, Kernel.Laplacian()))

function image_shot(idx, d_obs, q, model, origin; fmin=5f0, fmax=40f0)
    # Get current shot and project geometry
    newshot, newq = get_subset(data, q, origin, idx, fmin, fmax)
    DG = judiTimeGain(newshot.geometry, 1) # t^2
    Ml = judiDataMute(newq.geometry, newshot.geometry)
    newshot = Ml*DG*newshot

    # Setup operators with projected objects
    opt = Options(subsampling_factor=6, IC=ic_mode(model), limit_m=true,
                  free_surface=fs, space_order=16)  # ~40 GB of memory per source without subsampling
    M = judiModeling(model, newq.geometry, newshot.geometry; options=opt)

    # Imaging
    J = judiJacobian(M, newq)
    I = judiIllumination(J; mode="uv")
    rtm = J'*newshot
    
    return rtm, I
end


f32_read(x) = convert(Matrix{Float32}, transpose(segy_read(x).data))

if casenum == "077"
    vp = f32_read("$(data_path)/galactic-$(casenum)-Models/galactic-$(casenum)-FINAL-Vp.sgy") ./ 1f3
    vs = f32_read("$(data_path)/galactic-$(casenum)-Models/galactic-$(casenum)-FINAL-Vs.sgy") ./ 1f3
    vs[vs .== 1.0f-13] .= 0f0
    den = pyconvert(Matrix{Float32}, vtr.VTRModel("$(data_path)/galactic-$(casenum)-Models/galactic-$(casenum)-FINAL-Density.vtr").arrays[0]) ./ 1f3
    model_kw = Dict(:rho=> den, :vs=>vs)
    densname = "InvDensity"
else
    vp = f32_read("$(data_path)/galactic-$(casenum)-Models/galactic-$(casenum)-FINAL-Vp.sgy") ./ 1f3
    eps = f32_read("$(data_path)/galactic-$(casenum)-Models/galactic-$(casenum)-TrueEpsilon.sgy")
    del = f32_read("$(data_path)/galactic-$(casenum)-Models/galactic-$(casenum)-TrueDelta.sgy")
    del[del .> eps] .= .99f0 .* eps[del .> eps]

    if gardner
        den = f32_read("$(data_path)/galactic-$(casenum)-Models/galactic-$(casenum)-FINAL-Density-with-DenChange.sgy") ./ 1f3
        densname = "InvDensity"
    else
        den = pyconvert(Matrix{Float32}, vtr.VTRModel("$(data_path)/galactic-$(casenum)-Models/galactic-$(casenum)-FINAL-Density.vtr").arrays[0]) ./ 1f3
        densname = "Gardner"
    end
    model_kw = Dict(:rho=> den, :epsilon=>eps, :delta=>del)
end

n = size(vp)
d = (12.5f0, 12.5f0)

figure(figsize=(30, 12));
plot_velocity(vp', d; d_scale=0, new_fig=false, name="VP $(casenum)", save="$(figpath)/vp")
plotclose("all")

# Actual origin from original data
originm, o = get_orig(vp, "$(data_path)W22GAL_LINE$(linenum)_FTPreSTM_MigVel.segy")
@show o, originm

model = Model(n, d, (o, 0f0), vp.^(-2f0); model_kw...)
params = Dict(JUDI._params(model))

Tm = judiTopmute(model; taperwidth=0)

rtm = PhysicalParameter(model.G)
Ilu = PhysicalParameter(model.G)
Ilv = PhysicalParameter(model.G)

filename = "$(data_path)/galactic-$(casenum)-Models/rtm_$(densname)_125_$(fmax)$(multiples).bin"
batch_size = 20

if isfile(filename)
    rtm_prev = deserialize(filename)
    rtm = PhysicalParameter(rtm_prev.rtm, rtm_prev.spacing, rtm_prev.origin)
    Ilu = PhysicalParameter(rtm_prev.Ilu, rtm_prev.spacing, rtm_prev.origin)
    Ilv = PhysicalParameter(rtm_prev.Ilv, rtm_prev.spacing, rtm_prev.origin)
    iter = rtm_prev.iter + batch_size
else
    iter = 1
end

# idxsrc = vcat(collect.(i:100:data.nsrc for i in 1:100)...)
Random.seed!(data.nsrc)
idxsrc = randperm(data.nsrc)
last_src = data.nsrc

for it=iter:batch_size:last_src
    ite = min(data.nsrc, it+batch_size-1)
    i = idxsrc[it:ite]
    flush(stdout)
    t1 = @elapsed begin
        rtm_loc, Iloc = image_shot(i, data, q, model, originm; fmin=5f0, fmax=Float32(fmax))
        rtm_loc = Tm * rtm_loc
        if any(isnan.(rtm_loc.data))
            @info "$i has NaN"
            break
        end
        rtm .+= rtm_loc
        Ilu .+= Iloc.illums["u"]
        Ilv .+= Iloc.illums["v"]
        serialize(filename, (rtm=rtm.data, spacing=rtm.d, origin=rtm.o,
                             Ilu=Ilu.data, Ilv=Ilv.data, iter=it, params...))
    end
    println("Shot $(it) of $(data.nsrc) in: $(trunc(t1; digits=3)) s")

    Il = Ilu.^(0.5f0) .* Ilv
    rtm_plot = post_process(rtm.data, model) ./ (Il.data .+ 1f-6)
    figure(figsize=(30, 12));
    plot_simage(rtm_plot', rtm.d; cmap="cet_CET_L1", d_scale=0, perc=98,  new_fig=false, name="RTM $(casenum) $(densname)", save="$(figpath)/rtm_$(casenum)_$(densname)$(multiples)")
    plotclose("all")

end
