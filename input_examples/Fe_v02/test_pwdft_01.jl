using PWDFT

include("funcs_read.jl")
include("funcs_pwdft.jl")

function prepare_Rhoe_jdftx(Ham; do_transpose=false)
    n1, n2 = load_eVars_n(; do_transpose=do_transpose)
    Npoints = size(n1, 1)
    if sum(abs.(n2)) < 10*eps()
        Nspin = 1
    else
        Nspin = 2
    end
    @info "Npoints = $(Npoints)"
    @assert Npoints == prod(Ham.pw.Ns)
    @assert Nspin == Ham.electrons.Nspin

    Rhoe = zeros(Float64, Npoints, Nspin)
    if Nspin == 2
        Rhoe[:,1] = n1[:]
        Rhoe[:,2] = n2[:]
    else
        @assert Nspin == 1
        Rhoe[:,1] = n1[:]
    end

    return Rhoe
end

#=
function main()
    Ham = prepare_Ham();
    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns);

    Rhoe_jdftx = prepare_Rhoe_jdftx(Ham);

    Ham.electrons.Focc[:,:] = prepare_Focc(Ham);
    
    Rhoe = similar(Ham.rhoe);
    psiks = prepare_psiks(Ham);
    calc_rhoe!(Ham, psiks, Rhoe);

    Rhoe_ordered = similar(Ham.rhoe);
    psiks_new = reorder_psiks(Ham, psiks);
    calc_rhoe!(Ham, psiks_new, Rhoe_ordered);

end
=#

