using Random
using PWDFT

function main()

    Random.seed!(1234)
    
    atoms = Atoms(xyz_string_frac=
        """
        1

        Fe  0.0  0.0  0.0
        """, in_bohr=true,
        LatVecs = gen_lattice_bcc(5.42) )

    pspfiles = ["../../build/pseudopotentials/SG15/Fe_ONCV_LDA.upf"]
    ecutwfc = 20.0

    Ham = Hamiltonian(
        atoms,
        pspfiles,
        ecutwfc,
        meshk=[3,3,3],
        extra_states=4,
        Nspin=2
    );
    
    KS_solve_SCF!(
        Ham,
        mix_method="anderson",
        betamix=0.1,
        starting_magnetization=[0.1],
        use_smearing=true
    )

end

main()
