# This files contains structure that modelize functional space on which can be defined Hamiltonian

######## L²(Z) ########

struct L2Z
    Nmin::Int
    Nmax::Int
end

L2N(N::Int) = L2Z(0, N)

@inline Nmin(l2Z::L2Z) = l2Z.Nmin
@inline Nmax(l2Z::L2Z) = l2Z.Nmax

_zero(l2Z::L2Z) = zero(lattice(l2Z))

######## Complex Space ########

struct ComplexSpace end

const C = ComplexSpace()

_zero(::ComplexSpace) = 0.0 + 0.0*im

######## Lattice Function ########

lattice(l2Z::L2Z) = Nmin(l2Z):Nmax(l2Z)
lattice(::ComplexSpace) = nothing
lattice(t::Tuple) = Tuple([lattice(e) for e in t])

######## Base Generation ########

base(l2Z::L2Z, n::Int) = [n == i ? 1 : 0 for i in lattice(l2Z)]
base(l2Z::L2Z) = [base(l2Z, n) for n in lattice(l2Z)]
base(::ComplexSpace)  = 1.0 + 0.0*im

function base(t::Tuple)
    z = [_zero(e) for e in t]
    b = [base(e) for e in t]
    result = Vector{typeof(z)}() #[[_zero(e) for e in t ] for i in 1:sum(length.(b))]#Vector{typeof(z)}()
    for i in eachindex(t)
        for eb in b[i]
            push!(result, [j == i ? eb : z[j] for j in eachindex(t)])
        end
    end
    result
end