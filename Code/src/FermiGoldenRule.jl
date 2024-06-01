module FermiGoldenRule

    using Plots             # For plotting results
    using LaTeXStrings      # For the utilisation of Latex in Plots
    using Interpolations    # For interpolation used in the animation
    using LinearAlgebra     # For matrix computation
    using SparseArrays      # For matrix computation with sparsity
    using Distances         # For distances
    
    include("utils.jl")

    export relative_error_norm, proba
    include("analyse.jl")

    export L2Z, L2N, C, lattice, base 
    include("structure.jl")

    export dynamics
    include("dynamics.jl")

    export create_matrix
    include("model.jl")

    export animation, plot_proba_l2Z
    include("plot.jl")

end