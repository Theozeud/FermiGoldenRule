using Plots             # For plotting results
using LaTeXStrings      # For the utilisation of Latex in Plots
using LinearAlgebra     # For matrix computation
using SparseArrays      # For matrix computation with sparsity

using FermiGoldenRule

# Modelisation of the Hamiltonian
A(N::Int) = SymTridiagonal(zeros(N), ones(N-1))

δ(R::Int, Nmin::Int, Nmax::Int) = [R==i ? 1 : 0 for i in Nmin:Nmax] #Nmin:Nmax

h(E::Real, ϵ::Real, R₀::Int, Nmax::Int, Nmin::Int) = (z::Number, ϕ::Vector) -> (E*z + ϵ * ϕ[-Nmin+1+R₀], ϵ*z*δ(R₀,Nmin, Nmax)+A(length(ϕ))*ϕ)

# Initial condition
z₀ = 1.0
ϕ₀fun = (z₀, (N::Int) -> 0)

# Choice of model parameters
const R₀ = 0
E = 0
ϵ = 0.1

# Parameters for the approximation
ArrayNmax = [25,50,100,200] # Upper bound of the lattice
Nₜ = 2500                    # Definition of the number of timestep
T = 100                     # Duration of the simulation
timeT = LinRange(0,T,Nₜ+1)   # Time Lattice

# Initial condition on the lattiche
Arrayϕ₀ = [[ϕ₀fun[1], ϕ₀fun[2].(-Nmax:Nmax )...] for Nmax in ArrayNmax]


# Fermi's Golden Rule for the 1 site connected 1D chain model depending on the initial energy eigen
# Computation of the dynamics for different E and fixed ϵ
println("Beginning of the Fermi's Golden rule for the 1 site connected 1D chain model")

ϵ = 0.5
E = 0
ArrayNmaxsol= [dynamics(h(E,ϵ,R₀,Nmax,-Nmax), ϕ₀fun, (C,L2Z(-Nmax,Nmax)); Nₜ = Nₜ, T = T) for Nmax in ArrayNmax]

println("Dynamics computed !")

probaNmax = [proba(Arrayϕ₀[i], ArrayNmaxsol[i]) for i in eachindex(ArrayNmaxsol)]
plt = plot(size = (900,600), margin = 1Plots.cm, legendfontsize=14, legend = :bottomright, titlefontsize=14,
guidefontsize=14, tickfontsize=14)
for (proba,Nmax) in zip(probaNmax, ArrayNmax)
    plot!(plt, timeT,proba, label = latexstring("N = ",Nmax),lw = 5) 
end
xlabel!(L"t")
ylabel!(L"|⟨ϕ_0,ϕ(t)⟩|^2")
title!("Règle d'or de Fermi selon la valeur de N pour la chaine 1D pour ϵ = "*string(ϵ))
savefig(plt,"image/Règle d'or de Fermi pour la chaine 1D pour ϵ = 0,5 pour différentes valeurs de N")

println("Plot saved !")