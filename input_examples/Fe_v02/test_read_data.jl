using LinearAlgebra
using Serialization

function load_psiks()

    # Hardcoded parameters
    Nkspin = 8
    Nstates = 12
    Npw = [321, 350, 333, 334, 321, 350, 333, 334]

    Ndata = sum(Nstates * Npw)
    raw_data = Vector{ComplexF64}(undef, Ndata)

    file = open("eVars_C.dat", "r")
    read!(file, raw_data);
    close(file)

    psiks = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    idx_start = 1
    for ikspin in 1:Nkspin
        Ndata_k = Npw[ikspin]*Nstates
        idx_stop = idx_start + Ndata_k - 1
        psiks[ikspin] = zeros(ComplexF64, Npw[ikspin], Nstates)
        psiks[ikspin][:,:] = reshape(
            raw_data[idx_start:idx_stop], Npw[ikspin], Nstates
        )
        idx_start += Ndata_k
    end

    return psiks
end


function load_Hsub()
    Nkspin = 8
    Nstates = 12

    Ndata = sum(Nstates*Nstates * Nkspin)
    raw_data = Vector{ComplexF64}(undef, Ndata)

    file = open("eVars_Hsub.dat", "r")
    read!(file, raw_data);
    close(file)

    Hsub = Vector{Matrix{ComplexF64}}(undef, Nkspin)
    idx_start = 1
    for ikspin in 1:Nkspin
        Ndata_k = Nstates*Nstates
        idx_stop = idx_start + Ndata_k - 1
        Hsub[ikspin] = zeros(ComplexF64, Nstates, Nstates)
        Hsub[ikspin][:,:] = reshape(
            raw_data[idx_start:idx_stop], Nstates, Nstates
        ) # need transpose
        idx_start += Ndata_k
    end
    return Hsub
end

using Serialization
function main()
    #
    psiks = load_psiks()
    serialize("psiks.jldat", psiks)
    #
    Hsub = load_Hsub()
    serialize("Hsub.jldat", psiks)
end


