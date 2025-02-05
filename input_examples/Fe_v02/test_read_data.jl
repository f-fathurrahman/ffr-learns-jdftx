
using Serialization
function main()
    #
    psiks = load_psiks()
    serialize("psiks.jldat", psiks)
    #
    Hsub = load_Hsub()
    serialize("Hsub.jldat", Hsub)
    #
    Hsub_after = load_Hsub(filename="eVars_Hsub_after.dat")
    serialize("Hsub_after.jldat", Hsub_after)
end

#main()
