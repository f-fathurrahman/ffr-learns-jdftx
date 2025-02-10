function get_harcoded_params()
    # Hardcoded parameters
    Nkspin = 8
    Nstates = 7
    Npw = [321, 350, 333, 334, 321, 350, 333, 334]
    Ns = (20, 20, 20)
    return (;
        :Nkspin => Nkspin,
        :Nstates => Nstates,
        :Npw => Npw,
        :Ns => Ns
    )
end