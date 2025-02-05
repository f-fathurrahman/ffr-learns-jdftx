using LinearAlgebra
using Serialization

function get_harcoded_params()
    # Hardcoded parameters
    Nkspin = 8
    Nstates = 12
    Npw = [321, 350, 333, 334, 321, 350, 333, 334]
    return Nkspin, Nstates, Npw
end


function load_psiks()

    Nkspin, Nstates, Npw = get_harcoded_params()

    Ndata = sum(Nstates * Npw)
    raw_data = Vector{ComplexF64}(undef, Ndata)

    file = open("eVars_C.bindat", "r")
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


function load_Hsub(; filename="eVars_Hsub.bindat")

    Nkspin, Nstates, _ = get_harcoded_params()

    Ndata = Nstates * Nstates * Nkspin
    raw_data = Vector{ComplexF64}(undef, Ndata)

    file = open(filename, "r")
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
        ) # need transpose?
        idx_start += Ndata_k
    end
    return Hsub
end


function load_diagMatrix_vector(filename, T)

    Nkspin, Nstates, _ = get_harcoded_params()

    Ndata = Nstates * Nkspin
    raw_data = Vector{T}(undef, Ndata)

    file = open(filename, "r")
    read!(file, raw_data);
    close(file)

    D = Vector{Vector{T}}(undef, Nkspin)
    idx_start = 1
    for ikspin in 1:Nkspin
        Ndata_k = Nstates
        idx_stop = idx_start + Ndata_k - 1
        D[ikspin] = zeros(T, Nstates)
        D[ikspin][:] = raw_data[idx_start:idx_stop]
        idx_start += Ndata_k
    end
    return D
end

function load_eVars_n()
    n1 = readdlm("eVars_n_1.dat")
    if isfile("eVars_n_2.dat")
        n2 = readdlm("eVars_n_2.dat")
    else
        n2 = zeros(eltype(n1), size(n1))
    end
    return n1, n2
end
