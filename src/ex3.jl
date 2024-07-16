include("utils.jl")
using Plots

f0 = 5e9
z0 = 50
zLoad = 20 - 30im
zOpen = Inf 

lengthsA = [0.11, 0.125, 0.185]
lengthsB = [0.372, 0.125, 0.461]

function reflectionCoefCircuit(freq::Real, lengths) 
    newLengths = newLength.(lengths, freq, f0)

    zFirst = zInput(z0, zOpen, newLengths[1])
    z1 = parallel(zLoad, zFirst)

    z2 = zInput(z0, z1, newLengths[2]) 
    zSecond = zInput(z0, zOpen, newLengths[3])
    zTotal = parallel(z2, zSecond) 

    gamma = abs(reflectionCoef(z0, zTotal))
    return gamma
end

acceptableSWR = 2
acceptableGamma = (acceptableSWR - 1)/(acceptableSWR + 1) 
N = 1000
freqSweep = 0:2*f0/N:2*f0

gammaLimit = acceptableGamma*ones(N+1)
# small stub lengths
gammaA = reflectionCoefCircuit.(freqSweep, Ref(lengthsA))
gammaB = reflectionCoefCircuit.(freqSweep, Ref(lengthsB))

reflectionPlot = plot(freqSweep ./ 1e9, gammaA, label="small stub lengths")
plot!(freqSweep ./ 1e9, gammaB, color="red", label="big stub lengths") 
plot!(freqSweep ./ 1e9, gammaLimit, linestyle=:dot, color="black", label="0.33 bound")

plot!(xlabel="Frequency [GHz]", ylabel ="|Γ|", legend=:bottomleft)

