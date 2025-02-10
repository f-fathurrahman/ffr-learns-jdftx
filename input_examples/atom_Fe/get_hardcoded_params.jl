function get_harcoded_params()
    # Hardcoded parameters
    Nkspin = 2
    Nstates = 7
    Npw = [4337, 4337]
    Ns = (45, 45, 45)
    return (;
        :Nkspin => Nkspin,
        :Nstates => Nstates,
        :Npw => Npw,
        :Ns => Ns
    )
end
