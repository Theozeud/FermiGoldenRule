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
E = 0.1
ϵ = 0.0

# Parameters for the approximation
Nmax = 65                   # Upper bound of the lattice
Nmin = - Nmax               # Lower bound of the lattice    
Nₜ = 250                     # Definition of the number of timestep
T = 25                      # Duration of the simulation
timeT = LinRange(0,T,Nₜ+1)   # Time Lattice

# Computation of the dynamics
sol = dynamics(h(E, ϵ, R₀, Nmax, Nmin), ϕ₀fun, (C,L2Z(Nmin,Nmax)); T = T, Nₜ = Nₜ)

# Animation
info = (Nmax = Nmax, Nmin = Nmin)
animation(plot_simu, sol; trace =  true, d = 5, info = info)

