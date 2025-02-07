

filename = "iGarr.bindat"
Ndata = 3*691
raw_data = Vector{Int32}(undef, Ndata)

file = open(filename, "r")
read!(file, raw_data)
close(file)

