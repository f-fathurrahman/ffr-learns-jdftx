using Printf
using PWDFT

include("funcs_read.jl")

Ham, pwinput = init_Ham_from_pwinput(filename="PWINPUT");

Nkpt = Ham.pw.gvecw.kpoints.Nkpt
ik = 1;
ispin = 1;
ikspin = ik + (ispin-1)*Nkpt;

RecVecs = Ham.pw.RecVecs;
b1 = RecVecs[:,1];
b2 = RecVecs[:,2];
b3 = RecVecs[:,3];
# G_vec = b1*gi + b2*gj + b3*gk
kvec = Ham.pw.gvecw.kpoints.k[:,ik];
Ngw = Ham.pw.gvecw.Ngw;
idx_gw2g_k = Ham.pw.gvecw.idx_gw2g[ik];
G = Ham.pw.gvec.G;

iGarr = load_iGarr_ikspin(ikspin);

ecutwfc = 20.0
Ek_list = zeros(Float64, Ngw[ik])
for igw_jdftx in 1:Ngw[ik]
    gi, gj, gk = iGarr[igw_jdftx]
    Gkvec = b1*gi + b2*gj + b3*gk .+ kvec
    Ek = 0.5*dot(Gkvec,Gkvec)
    Ek_list[igw_jdftx] = Ek
    @printf("miller_jdftx=(%4d,%4d,%4d) Ek=%18.10f\n", gi, gj, gk, Ek)
    if Ek > ecutwfc
        @warn "Larger than ecutwfc"
    end
end


