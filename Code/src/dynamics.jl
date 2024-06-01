# This files contains methods to run the dynamics of a system

init(ϕ₀fun::Tuple, Lattice)  = vcat([apply(ϕ₀fun[i],Lattice[i]) for i in eachindex(Lattice)]...)
init(ϕ₀fun::Base.Callable, Lattice)  = apply(ϕ₀fun,Lattice)
 
function dynamics(Hfun::Base.Callable, ϕ₀fun::Tuple, Infos_discr::Tuple; T::Real = 4, Nₜ::Int = 1000)

    # INPUTS
    #   - Hfun    : Hamiltonian of the model,
    #   - ϕ₀fun   : Initial condition at time 0,
    #   - T       : Duration defining the temporal domain [0,T] on wich we perform the simulation,
    #   - Nₜ       : Number of points used for an uniform discretization of the time.

    # OUTPUTS
    #   - matrix of dimension Nₓ×(Nₜ+1) the columns of which correspond to the solution at each time beginning at the initial state.
   
    # lattice
    Lattice = lattice(Infos_discr)

    # Initial Condition
    ϕ₀ = init(ϕ₀fun, Lattice) 

    # Matrix version of the approximate Hamiltonian
    H = create_matrix(Hfun, Infos_discr)

    dynamics(H, ϕ₀; T = T, Nₜ = Nₜ)
end

function dynamics(H::Matrix, ϕ₀::Vector; T::Real = 4, Nₜ::Int = 1000)

    # INPUTS
    #   - H        : Hamiltonian of the model,
    #   - ϕ₀       : Initial condition at time 0,
    #   - T        : Duration defining the temporal domain [0,T] on wich we perform the simulation,
    #   - Nₜ        : Number of points used for an uniform discretization of the time.

    # OUTPUTS
    #   - matrix of dimension Nₓ×(Nₜ+1) the columns of which correspond to the solution at each time beginning at the initial state.
   
    # Timestep
    Δt = eltype(ϕ₀)(T/Nₜ)

    # Matrix to reach next each step
    stepH = exp(-im * Δt * H)

    # Matrix to store the solution
    sol = zeros(Complex, length(ϕ₀),Nₜ+1)
    sol[:,1] = ϕ₀
    
    # Loop to compute the solution for each time step
    for i in 2:Nₜ+1
        @inbounds sol[:,i] = stepH * sol[:,i-1]
    end
    sol
end
