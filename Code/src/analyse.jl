# Some Functionality to analyse the systems

# Relative error made on the norm of the solution
relative_error_norm(sol::Matrix{Complex}, i::Int) =  abs(norm(sol[:,i])/norm(sol[:,1])-1)*100
relative_error_norm(Sol::AbstractVector, i::Int) =  [relative_error_norm(sol,i) for sol in Sol]
relative_error_norm(sol::Matrix{Complex}, I::AbstractVector = axes(sol)[2]) = max([relative_error_norm(sol, i) for i in I]...)
relative_error_norm(Sol::AbstractVector, I::AbstractVector) = [relative_error_norm(sol,I) for sol in Sol]
relative_error_norm(Sol::AbstractVector) = [relative_error_norm(sol) for sol in Sol]

# Functions to compute the probability of presence
proba(ϕ::Vector, ϕ₀::Vector) = abs(dot(ϕ, ϕ₀))^2
proba(ϕ₀::Vector, ϕ::Matrix) = vcat([proba(ϕ[:,i], ϕ₀) for i in axes(ϕ)[2]]...)

