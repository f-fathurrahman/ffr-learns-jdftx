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
