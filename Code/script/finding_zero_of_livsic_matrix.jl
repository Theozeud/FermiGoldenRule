using Plots             # For plotting results
using LaTeXStrings      # For the utilisation of Latex in Plots
using LinearAlgebra     # For matrix computation
using Optim             # For optimization function
using Distances         # For Distance function

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
E = 0.0
ϵ = 0.1

# Parameters for the approximation
Nmax = 100                  # Upper bound of the lattice
Nmin = - Nmax               # Lower bound of the lattice    
Lattice = Nmin:Nmax         # Definition of the discrete lattice
Nₜ = 250                     # Definition of the number of timestep
T = 25                      # Duration of the simulation
timeT = LinRange(0,T,Nₜ+1)   # Time Lattice

# Initial condition on the lattiche
ϕ₀ = [ϕ₀fun[1], ϕ₀fun[2].(Lattice)...]

# Computation of the dynamics
sol = dynamics(h(E, ϵ, R₀, Nmax, Nmin), ϕ₀fun, (C,L2Z(Nmin,Nmax)); T = T, Nₜ = Nₜ)

# To find the coefficient Γ(ϵ), we minimze the function equal to the solution minus exp(-2Γ(ϵ)tϵ^2)

# Loss function
loss(Γ, sol, time) = norm(sqeuclidean(proba(ϕ₀,sol),exp.(-2*Γ[1].*time*ϵ^2)))^2

# Init guess for the optimization
E_Array = LinRange(-2,2,100)
Γ_ϵ = []
Γ₀ = 1.0

@time for e in E_Array
    sol_E = dynamics(h(e, ϵ, R₀, Nmax, Nmin), ϕ₀fun, (C,L2Z(Nmin,Nmax)); T = T, Nₜ = Nₜ)
    opt = optimize(Γ -> loss(Γ, sol_E, timeT), [Γ₀]; iterations=100)
    res = Optim.minimizer(opt)
    push!(Γ_ϵ, res[1])
end

plt = plot(E_Array, 1 ./ (2 .* sqrt.(1 .- E_Array .^2 ./ 4)), label = "Théorique", lw=5, ls = :dot)
plot!(E_Array, Γ_ϵ, size = (900,600), margin = 1Plots.cm, legendfontsize=14, legend = :bottomright, titlefontsize=14,
guidefontsize=14, tickfontsize=14, lw=5, label = "Numérique")
xlabel!(L"E")
ylabel!(L"-Γ(ϵ)/ϵ^2")
title!("-Γ(ϵ)/ϵ^2 selon la valeur de la valeur propre initial E")
savefig(plt,"image/Γ(ϵ) selon la valeur de la valeur propre initial E")