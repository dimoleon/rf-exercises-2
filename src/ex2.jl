include("radiation.jl")
using PhysicalConstants.CODATA2018: c_0, ε_0, μ_0

plotlyjs()

I = 1
f = 1e9
λ = c_0.val/f
k = 2π/λ
r = 1
N = 8

dlist = [0.25λ, 0.5λ, 0.75λ]

function E0(θ) 
    res = abs(im*60*I*cis(-k*r)*cosd(90*cosd(θ))/(r*sind(θ)))
    if isnan(res)
        res = 0
    end
    return res
end

function E(ϕ, θ, d)
    res = 2*abs(E0(θ))*abs(sum(cos(m*k*d*cosd(ϕ)*sind(θ)) for m ∈ 0.5:(0.5(N-1))))
    return res
end

function Dmultidirectional(d) 
    return 2*N*d/λ
end

# change index on the next line for the other diagrams of the exercise
curr_d = dlist[3]
Efix(ϕ, θ) = E(ϕ, θ, curr_d)
maxSteps = 361

phi = range(0, stop=360, length=maxSteps)
p = plot_xy(Efix, phi)

theta = range(0, stop=360, length=maxSteps)
p = plot_xz(Efix, theta)

phi = range(0, stop=360, length=maxSteps)
theta = range(0, stop=180, length=maxSteps)
p = plot_3d(Efix, phi, theta)

dir = directivity(Efix)
dir_mult = Dmultidirectional(curr_d)
println("directivity: $dir")
println("multi directonial: $dir_mult")