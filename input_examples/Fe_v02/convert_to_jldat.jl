using Serialization
using PWDFT

include("funcs_read.jl")
include("funcs_pwdft.jl")

function main_v1()
    #
    psiks = load_psiks()
    serialize("psiks.jldat", psiks)
    #
    Hsub = load_Hsub()
    serialize("Hsub.jldat", Hsub)
    #
    Hsub_after = load_Hsub(filename="eVars_Hsub_after.bindat")
    serialize("Hsub_after.jldat", Hsub_after)
end
#main_v1()


# Serialize psiks (reordered) and Focc
function main()
    Ham = prepare_Ham()
    Focc = prepare_Focc(Ham)
    psiks = prepare_psiks(Ham);
    psiks_new = reorder_psiks(Ham, psiks);
    serialize("psiks.jldat", psiks_new)
    serialize("Focc.jldat", Focc)
    #
    Hsub = load_Hsub()
    serialize("Hsub.jldat", Hsub)
    #
    g_Haux = load_matrices_Nkspin_Nstates("g_Haux.bindat")
    serialize("g_Haux.jldat", g_Haux)
    Kg_Haux = load_matrices_Nkspin_Nstates("Kg_Haux.bindat")
    serialize("Kg_Haux.jldat", Kg_Haux)
end
main()


