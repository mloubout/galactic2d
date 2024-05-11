using SegyIO, JUDI, SlimPlotting, Interpolations, Serialization, SourceEstimation, LinearAlgebra, JLD2
close("all")
data_path = "/data/galactic2D/"
linenum = 18
include("utils.jl")

data, q = load_slice(linenum, "/data/galactic2D/GAL_FullArray_FFSig_Ver2_68pt5ms_0pt5ms.sgy")

function image_shot(idx, d_obs, q, model, origin)
    # Get current shot and project geometry
    newshot = get_data(d_obs[idx]; rel_origin=origin, project="2d")
    newq = get_data(q[idx]; rel_origin=origin, project="2d")
    Fq = judiFilter(newq, 3f0, 30f0)
    newq = Fq * newq
    # Setup operators with projected objects
    opt = Options(subsampling_factor=6, IC="isic", limit_m=true, free_surface=true, space_order=16)  # ~40 GB of memory per source without subsampling
    M = judiModeling(model, newq.geometry, newshot.geometry; options=opt)
    # Inversion
    # F = judiFilter(newshot, 3f0, 50f0)
    J = judiJacobian(M, newq)
    # d_syn = M*newq
    # P, qest = estimate_source!(d_syn, F*newshot, newq; taperwidth=300, ntrace=50,  nt=300, lambda=1f-1, beta=1f-6)

    # r = P'*(d_syn - newshot)
    # J = judiJacobian(M, qest)
    I = judiIllumination(J; mode="v")
    rtm = J'*newshot
    
    return rtm, I.illums["v"]
end


start_vp = "$(data_path)galactic-A85-StartVp.sgy"
xwi_vp = "$(data_path)galactic-A85-CP00054-Vp.sgy"

for (vp, name) in [(start_vp, "xwi_start"), (xwi_vp, "xwi_54")]

    model, originm  = read_model("$(data_path)W22GAL_LINE$(linenum)_FTPreSTM_MigVel.segy"; d=12.5,
                                vals=vp)

    rtm = similar(model.m)
    fill!(rtm, 0)
    Il = similar(model.m)
    fill!(Il, 0.)

    filename = "$(data_path)rtm_line_18_$(name)_125.bin"

    if isfile(filename)
        rtm, _, dloc, oloc, Il, iter = deserialize(filename)
        rtm = PhysicalParameter(rtm, dloc, oloc)
        Il = PhysicalParameter(Il, dloc, oloc)
        iter += 1
    else
        iter = 1
    end

    for i=iter:data.nsrc
        flush(stdout)
        t1 = @elapsed begin
            rtm_loc, Il_loc = image_shot(i, data, q, model, originm)
            rtm .+= rtm_loc
            Il .+= Il_loc
            serialize(filename, (rtm=rtm.data, m=model.m.data, spacing=rtm.d, origin=rtm.o, Il=Il.data, iter=i))
        end
        println("Shot $(i) of $(data.nsrc) in: $(trunc(t1; digits=3)) s")
    end

end
