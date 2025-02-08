using LinearAlgebra
using Serialization
using DelimitedFiles

function get_harcoded_params()
    # Hardcoded parameters
    Nkspin = 8
    Nstates = 11
    Npw = [691, 683, 678, 675, 691, 683, 678, 675]
    Ns = (24, 24, 24)
    return (;
        :Nkspin => Nkspin,
        :Nstates => Nstates,
        :Npw => Npw,
        :Ns => Ns
    )
end

function load_iGarr_ikspin(ikspin)
    (; Npw) = get_harcoded_params()

    filename = "iGarr_" * string(ikspin) * ".bindat"
    Ndata = 3*Npw[ikspin]
    raw_data = Vector{Int32}(undef, Ndata)
    #
    file = open(filename, "r")
    read!(file, raw_data)
    close(file)
    #
    idx_miller = Vector{Tuple{Int64,Int64,Int64}}(undef, Npw[ikspin])
    istart = 1
    for i in 1:Npw[ikspin]
        idx_miller[i] = (raw_data[istart], raw_data[istart+1], raw_data[istart+2])
        istart += 3
    end
    return idx_miller
end

function load_psiks()

    (; Nkspin, Nstates, Npw) = get_harcoded_params()

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

    (; Nkspin, Nstates) = get_harcoded_params()

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

    (; Nkspin, Nstates) = get_harcoded_params()

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

function load_eVars_n(; do_transpose=false)

    (; Ns) = get_harcoded_params()

    n1 = readdlm("eVars_n_1.dat")
    @assert size(n1) == (prod(Ns),1)
    if isfile("eVars_n_2.dat")
        n2 = readdlm("eVars_n_2.dat")
        @assert size(n2) == (prod(Ns),1)
    else
        n2 = zeros(eltype(n1), size(n1))
    end

    if do_transpose
        n1 = permutedims(
            reshape(n1, Ns),
            (3,2,1)
        )
        n1 = reshape(n1, prod(Ns), 1)

        n2 = permutedims(
            reshape(n2, Ns),
            (3,2,1)
        )
        n2 = reshape(n2, prod(Ns), 1)
    end

    return n1, n2
end
