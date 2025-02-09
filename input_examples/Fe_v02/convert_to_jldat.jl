using Serialization

include("funcs_read.jl")

function main()
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

main()
