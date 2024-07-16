using PhysicalConstants.CODATA2018: c_0
using NLsolve

c0 = c_0.val
f0 = 2.5e9
λ0 = c0/f0
k0 = 2π/λ0
z0 = 50
d = 1.6e-3
er = 4.4
Q = 25

A = z0/60*sqrt(0.5(er + 1)) + (er-1)/(er +1)*(0.23 + 0.11/er)
w = 8exp(A)/(exp(2A) - 2)*d

eff = 0.5(er + 1) + 0.5(er - 1)/sqrt(1 + 12d/w)

λ = λ0/sqrt(eff)

l = 0.5λ

β = 2π*f0*sqrt(eff)/c0
α = β/(2Q)

function f!(F, x)
    F[1] = tan(2π*x[1]*sqrt(eff)*l/c0) + sqrt(α*c0/(2x[1]*sqrt(eff)))
end

res = nlsolve(f!, [f0])
f_resonance = res.zero[1]

β_res = β*f_resonance/f0
b_gap = -tan(β_res*l)

C_gap = b_gap/(2π*f_resonance*z0)


# (c) 
tandelta = α*2*sqrt(eff)*(er - 1)/(k0*er*(eff - 1))

lamda(f) = c0/(sqrt(eff)*f)
k(f) = 2π*f/c0
alpha(f) = k(f)*er*(eff - 1)*tandelta/(2*sqrt(eff)*(er -1))
alpha(f0)

function g!(F, y)
    F[1] = tan(π*f0*lamda(y[1])*sqrt(eff)/c0) + sqrt(complex(alpha(y[1])*c0/(2*f0*sqrt(eff))))
end

res_new = nlsolve(g!, [f0])
f_simple = res_new.zero[1]
l_simple = lamda(f_simple)/2

b_s = -tan(β*l_simple)
C_s = b_s/(2π*f_simple*z0)

