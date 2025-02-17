function prepare_Ham()
    Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT")
    return Ham
end

function get_idx_pwdft_G(Ham, ikspin)
    iGarr = load_iGarr_ikspin(ikspin)
    idx_g2miller = Ham.pw.gvec.idx_g2miller
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    if Nspin == 2
        if ikspin > Nkpt
            ik = ikspin - Nkpt
        else
            ik = ikspin
        end
    else
        @assert Nspin == 1
        ik = ikspin
    end
    Ngwk = Ham.pw.gvecw.Ngw[ik]
    idx_gw2g = zeros(Int64, Ngwk)
    for igw in 1:Ngwk
        ig1 = iGarr[igw]
        ig2 = findfirst(isequal(ig1), idx_g2miller)
        idx_gw2g[igw] = ig2
    end
    return idx_gw2g
end

function reorder_psiks(Ham, psiks)
    #
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Ngw = Ham.pw.gvecw.Ngw
    psiks_new = zeros_BlochWavefunc(Ham)
    #
    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        idx_gw2g_k = Ham.pw.gvecw.idx_gw2g[ik]
        idx_G = get_idx_pwdft_G(Ham, ikspin)
        for igw in 1:Ngw[ik]
            igw_pwdft = findfirst(isequal(idx_G[igw]), idx_gw2g_k);
            #@info "igw=$(igw) idx_G=$(idx_G[igw]) igw_pwdft=$(igw_pwdft)"
            psiks_new[ikspin][igw_pwdft,:] = psiks[ikspin][igw,:]
        end
    end
    return psiks_new
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