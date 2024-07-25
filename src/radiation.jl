using Plots

plotlyjs()
function plot_xy(field, phisweep)
    Ehor = field.(phisweep, 90)

    p = plot(deg2rad.(phisweep), Ehor, proj=:polar, title="Horizontal xy plane", legend=false)

    return p 
end

function plot_xz(field, thetasweep)
    Ever = field.(0, thetasweep)

    p = plot(deg2rad.(thetasweep), Ever, proj=:polar, title="Vertical xz plane", legend=false, gridlinewidth=2.5, guidefontsize=5, tickfontsize=5, ratio=1)

    return p 
end

function plot_3d(field, phisweep, thetasweep) 
    steps = length(thetasweep)
    s = field.(phisweep, permutedims(thetasweep))
    x = s .* (cosd.(phisweep)* sind.(thetasweep)')
    y = s .* (sind.(phisweep) * sind.(thetasweep)')
    z = s .* (ones(steps) * cosd.(thetasweep)')

    p = plot(surface(x, y, z, color=:diverging, surfacecolor=@.sqrt(x^2 + y^2 + z^2)), legend=false, xlabel="X", ylabel="Y", zlabel="Z", title="3D Radiation Diagram")

    maxs = maximum(s)
    xlims!(-maxs*1.2, maxs*1.2)
    ylims!(-maxs*1.2, maxs*1.2)
    zlims!(-maxs*1.2, maxs*1.2)

    return p 
end

function directivity(field)
    phi = range(0, stop=359, length=360)
    theta = range(0, stop=180, length=181)

    s = field.(phi, permutedims(theta))

    P = s.^2 
    Pmax = maximum(P) 
    
    Wr = sum(P * sind.(theta))
    dtheta = π/180
    dphi = 2π/360
    Wr = Wr * dtheta * dphi

    D = 4π * Pmax / Wr
    return D
end
