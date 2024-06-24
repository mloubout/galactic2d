using SegyIO, JUDI, SlimPlotting, Interpolations, Serialization, LinearAlgebra, JLD2, Random

fmax = 40

close("all")
data_path = "/data/galactic2D/"
figpath = "./images/rtm-$(fmax)-2405"
mkdir(figpath)
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
    
    return rtm, I.illums["v"], I.illums["u"]
end


start_vp = "$(data_path)galactic-StartVp-0524.sgy"
xwi_vp = "$(data_path)galactic-CP00054-Vp-0524.sgy"

for (vp, name) in [(xwi_vp, "xwi_54")] #, (start_vp, "xwi_start")]

    model, originm  = read_model("$(data_path)W22GAL_LINE$(linenum)_FTPreSTM_MigVel.segy"; d=12.5,
                                 vals=vp, density=true)

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

    filename = "$(data_path)rtm_line_18_0524_$(name)_125_$(fmax).bin"

    if isfile(filename)
        rtm, _, dloc, oloc, Ilu, Ilv, iter = deserialize(filename)
        rtm = PhysicalParameter(rtm, dloc, oloc)
        Ilu = PhysicalParameter(Ilu, dloc, oloc)
        Ilv = PhysicalParameter(Ilv, dloc, oloc)
        iter += 1
    else
        iter = 1
    end

    # idxsrc = vcat(collect.(i:100:data.nsrc for i in 1:100)...)
    Random.seed!(data.nsrc)
    idxsrc = randperm(data.nsrc)

    for it=iter:data.nsrc
        i = idxsrc[it]
        flush(stdout)
        t1 = @elapsed begin
            rtm_loc, Ilu_loc, Ilv_loc = image_shot(i, data, q, model, originm; fmin=5f0, fmax=Float32(fmax))
            rtm_loc = Tm * rtm_loc
            if any(isnan.(rtm_loc.data))
                @info "$i has NaN"
                continue
            end
            rtm .+= rtm_loc
            Ilu .+= Ilu_loc
            Ilv .+= Ilv_loc
            serialize(filename, (rtm=rtm.data, m=model.m.data, spacing=rtm.d, origin=rtm.o,
                                 Ilu=Ilu.data, Ilv=Ilv.data, iter=it))
        end
        println("Shot $(it) of $(data.nsrc) in: $(trunc(t1; digits=3)) s")
    end

end
