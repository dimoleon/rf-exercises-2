include("radiation.jl")
using PhysicalConstants.CODATA2018: c_0, ε_0, μ_0

f = 1e9
λ = c_0.val / f 
k = 2π/λ
r = 1

d1 = λ/sqrt(8)
d2 = 3λ/sqrt(8)

heights = [[d1, d1], [d1, d2], [d2, d2]] 

function E0(θ) 
    res = abs(im*60*I*cis(-k*r)*cosd(90*cosd(θ))/(r*sind(θ)))
    if isnan(res)
        res = 0
    end
    return res
end

function E(ϕ, θ, h) 
    papa = abs(sin(k*h[1]*cosd(ϕ)*sind(θ)) * sin(k*h[2]*sind(ϕ)*sind(θ)))
    return papa * 4 * E0(θ)
end

currh = heights[3]
Ecurr(ϕ, θ) = E(ϕ, θ, currh)
maxSteps = 181
phi = range(0, 90, length = maxSteps)
theta = range(0, 180, length = maxSteps)
p = plot_xy(Ecurr, phi)

p = plot_3d(Ecurr, phi, theta)