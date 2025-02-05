using PWDFT

include("funcs_read.jl")

function prepare_Ham()
    Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT")
    return Ham
end

function prepare_psiks(Ham)
    
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates
    CellVolume = Ham.pw.CellVolume

    psiks = load_psiks()
    Nkspin = length(psiks)
    @assert Nkspin == Nkpt*Nspin

    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        @assert Ham.pw.gvecw.Ngw[ik] == size(psiks[ikspin], 1)
        @assert Nstates == size(psiks[1], 2)
        #
        # Renormalize, assuming NCPP, no overlap operator
        psiks[ikspin][:,:] .*= sqrt(CellVolume)
    end
    return psiks
end


function prepare_Focc(Ham)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates

    Focc = load_diagMatrix_vector("eVars_F.bindat", Float64)
    Nkspin = length(Focc)

    @assert Nkpt*Nspin == Nkspin

    for ikspin in 1:Nkspin
        @assert length(Focc[ikspin]) == Nstates
    end

    # convert to matrix
    Focc_r = zeros(Float64, Nstates, Nkspin)
    for ikspin in 1:Nkspin
        Focc_r[:,ikspin] = Focc[ikspin]
    end

    return Focc_r
end

function prepare_Rhoe_jdftx(Ham)
    n1_, n2_ = load_eVars_n()
    Npoints = size(n1_, 1)
    if sum(abs.(n2_)) < 10*eps()
        Nspin = 1
    else
        Nspin = 2
    end
    @assert Npoints == prod(Ham.pw.Ns)
    @assert Nspin == Ham.electrons.Nspin

    n1 = permutedims(
        reshape(n1_, (20,20,20)),
        (3,2,1)
    )

    n2 = permutedims(
        reshape(n2_, (20,20,20)),
        (3,2,1)
    )

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
    Ham = prepare_Ham()
    psiks = prepare_psiks(Ham)
end
=#

