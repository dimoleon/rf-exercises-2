using PhysicalConstants.CODATA2018: c_0, ε_0, μ_0
using Printf
using NLsolve

a = 2.286e-2
b = 1.016e-2
d = 1.5e-3
L = 6e-2
f = 10e9

c0 = c_0.val
ε0 = ε_0.val
μ0 = μ_0.val
η0 = sqrt(μ0/ε0)
k0 = 2π*f/c0

fcut = c0/(2a)
bigroot = sqrt(1 - (fcut/f)^2) 

β0 = 2π*f*bigroot/c0

Z0 = η0/bigroot

Zin1 = 4.9678+43.9439im
Zin2 = 108.5347+202.0158im

ZA = Z0*(Zin1 - Z0*tan(β0*(L-d))im)/(Z0 - Zin1*tan(β0*(L - d))im)
ZB = Z0*(Zin2 - Z0*tan(β0*(L-2d))im)/(Z0 - Zin2*tan(β0*(L - 2d))im)

n = 0
γd = atanh(sqrt(2ZA/ZB - 1)) + π*n*im
γ = γd / d

α = real(γ)
β1 = imag(γ)

ϵr = (β1*c0/(2π*f))^2 + (fcut/f)^2

k1 = k0*sqrt(ϵr)

tanδ = 2α*β1/k1^2

@printf("1.a ϵ_r = %f\n", ϵr)
@printf("1.a tanδ = %f\n", tanδ)

# Now we only have the first measurement available

C1 = ZA/(2π*f*μ0*im*d)
C2 = ZB/(2π*f*μ0*im*2d)

function f!(F, x)
    z = x[1] + 1im * x[2]
    F[1] = real(tanh(z)/z - C1)
    F[2] = imag(tanh(z)/z - C1)
end

function g!(F, x)
    z = x[1] + 1im * x[2]
    F[1] = real(tanh(z)/z - C2)
    F[2] = imag(tanh(z)/z - C2)
end

res1 = nlsolve(f!, [real(γd), imag(γd)])
gd1 = complex(res1.zero[1], res1.zero[2])

res2 = nlsolve(g!, [0.1, 0.1])
gd2 = complex(res2.zero[1], res2.zero[2])/2

ab, bb = real(gd1/d), imag(gd1/d)
ϵrb = (bb*c0/(2π*f))^2 + (fcut/f)^2
k1 = k0*sqrt(ϵrb)

tanδ = 2ab*bb/k1^2
println(γ)
println(gd1/d)
println(gd2/d)
