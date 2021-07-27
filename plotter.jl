using XMLDict
using Plots
using LinearAlgebra
using LaTeXStrings
using Plots.PlotMeasures

function compute_kx(k)
    n_kpts = size(k, 2)
    # Allocate Arrays
    kx = zeros(n_kpts)
    # When symm[n] is true, it means k[n] is a high symm point
    symm = fill(false, n_kpts)
    symm[1] = true # First point is always high symmetry

    # Calculate stuff for the first step
    k_tmp = k[:,2] - k[:,1]
    dk_tmp = sqrt(k_tmp ⋅ k_tmp)

    # Temporary Arrays
    k1 = Vector{Float64}(undef, 3)
    k2 = Vector{Float64}(undef, 3)

    for n in 2:n_kpts
        k_tmp = k[:,n] - k[:,n-1] 
        dk = sqrt(k_tmp ⋅ k_tmp)
        # If there's a large jump, the path is discontinuous and we overlay the points
        if dk > 5*dk_tmp
            kx[n] = kx[n-1]
        elseif dk > 1e-5 # normal case: points separated by small delta
            kx[n] = kx[n-1] + dk
            dk_tmp = dk
        else # points separated by delta almost zero, add to path but don't save delta
            kx[n] = kx[n-1] + dk
        end
        # Then we determine symmetry
        if k[:,n] ⋅ k[:,n] < 1e-9 # gamma point is high symm
            symm[n] = true
        elseif n == n_kpts
            symm[n] = true # Last point is also a high symmetry point
        else
            k1 = k[:,n] - k[:,n-1]
            k2 = k[:,n+1] - k[:,n]
            cos_theta = (k1⋅k2)/sqrt(k1⋅k1)/sqrt(k2⋅k2)
            # If the path changes direction then it's a high symm point
            if abs(cos_theta - 1) > 1.0e-4
                symm[n] = true
            end
        end
    end
    
    symm_kx = kx[symm]
    symm_k = k[:,symm]

    for i in eachindex(symm_kx)
        @info "High Symmetry Point $(symm_k[:,i]) at kx = $(symm_kx[i])"
    end

    return kx, symm_kx
    
end

function parse_xml(path)
    xml_string = open(path) do file
        read(file, String)
    end

    # Read as XML dictionary
    d = xml_dict(xml_string)

    # These are the k-points in an array where each element is a string of 3 numbers
    kpoints = d["modeling"][""][10]["kpoints"]["varray"][1]["v"]
    nk = length(kpoints) # number of k-points
    k = Array{Float64}(undef, 3, nk)
    # Loop over k points and save them into k, in proper float format
    for n in 1:nk
        kpoint = kpoints[n]
        k[:,n] =  map(x->parse(Float64, x),split(kpoint))
    end

    # Access the eigenvalues (and occupations, lol)
    eigenvalues = d["modeling"][""][18]["calculation"]["eigenvalues"]["array"]["set"]["set"]["set"]
    nb = length(eigenvalues[1]["r"]) # Number of bands
    E = Array{Float64}(undef, nb, nk)
    occ = Array{Float64}(undef, nb, nk)
    for n in 1:nk
        for m in 1:nb
            eigenvalue = eigenvalues[n]["r"][m]
            E[m, n] = parse(Float64, split(eigenvalue)[1])
            occ[m, n] = parse(Float64, split(eigenvalue)[2])
        end
    end
    return k, E, occ
end

# This is not particularly accurate with how folks smear in VASP.
function compute_band_props(E, occ)
    # Occupied bands are those with 1 in the occupations
    occupied = findall(==(1), occ)
    occupied_bands = map(i->i[1], occupied)
    # Unoccupied bands are those with 0 in the occupations
    unoccupied = findall(==(0), occ)
    unoccupied_bands = map(i->i[1], unoccupied)
    # Maximum occupied and minimum unoccupied bands
    max_vb = maximum(occupied_bands)
    min_cb = minimum(unoccupied_bands)
    # Compute vbm and cbm and band gap
    vbm = maximum(E[max_vb, :])
    cbm = minimum(E[min_cb, :]) 
    Eg = (cbm - vbm)*1000
    @info "Band Gap: $Eg meV"
    @info "CBM: $cbm eV"
    @info "VBM: $vbm eV"

    return vbm, cbm, max_vb, min_cb
end

function plot_bands(path)
    (k, E, occ) = parse_xml(path)
    # Compute linearized k and plot
    kpath = [L"R", L"\Gamma", L"X", L"M", L"\Gamma"]
    (kx, symm)  = compute_kx(k)
    (vbm, cbm, max_vb, min_cb) = compute_band_props(E, occ)
    p = plot(kx, E[1:max_vb, :]'.- vbm, lc = :blue, lw=2, label="", fmt=:png, dpi=200)
    plot!(p, kx, E[min_cb:end, :]'.- vbm, lc = :orange, lw=2, ylim=(-6,6), label="")
    # Plot reference energy and symmetry lines
    plot!(p, [0], seriestype=:hline, linestyle=:dash, lw=1.5, label="", linecolor=:black)
    plot!(p, symm, seriestype=:vline, linestyle=:dashdot, lw=1, label="", linecolor=:black)
    plot!(p, ylabel = L"E - E_{\textrm{vbm}}", xticks = (symm, kpath), grid = false, xlim=(minimum(symm), maximum(symm)))
    plot!(p, xtickfontsize=16, ytickfontsize=14, yguidefontsize=16, margin=3mm)
    return p
end

p = plot_bands("vasprun.xml")