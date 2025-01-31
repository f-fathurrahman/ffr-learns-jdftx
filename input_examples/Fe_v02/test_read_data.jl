using LinearAlgebra

function test_main()

    Nkspin = 8
    Nstates = 11
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
        println("idx_start = ", idx_start)
        println("idx_stop = ", idx_stop)
        psiks[ikspin] = zeros(ComplexF64, Npw[ikspin], Nstates)
        psiks[ikspin][:,:] = reshape(
            raw_data[idx_start:idx_stop], Npw[ikspin], Nstates
        )
        idx_start += Ndata_k
    end

    return psiks
end

test_main()
