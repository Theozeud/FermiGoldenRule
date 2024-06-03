function animation(funplot::Base.Callable, ψ; d::Real = 5, fps::Int = 30, trace::Bool = false, info::NamedTuple)

    # INPUTS
    #   - funplot   : Function to plot,
    #   - ψ         : Solution Matrix given by the dynamics function,
    #   - d         : Duration wanted for the animation,
    #   - fps       : Number of frame per secondes,
    #   - trace     : Boolean to plot trace or not,
    #   - info      : Set of information to give to plot functions.

    # OUTPUTS
    #   - An animation.
    
    # Some assertions to check validity of parameters
    @assert d > 0
    @assert fps > 0
    
    # Total frame
    nbframe = fps * d
  
    # Time vector
    time = range(0,1,size(ψ)[2])
    
    # Storage of trace
    storage = Complex{Float64}[]
  
    anim = @animate for n in 1:nbframe

        # Actual frame
        t = n / (fps*d)
        # Interpolation
        interp_module = [linear_interpolation(time, abs.(ψ[i,:]))(t) for i in eachindex(ψ[:,1])]
        interp_arg    = [linear_interpolation(time, angle.(ψ[i,:]))(t) for i in eachindex(ψ[:,1])]
        interp_spher  = interp_module .* exp.(im.*interp_arg)

        Trace = funplot(interp_spher; info = info)

        if trace
            push!(storage, Trace)
            plot!(storage, color = :blue)
        end
        
    end
  
    gif(anim, fps = fps)
  end


  function plot_compoC(Φ::Number, t::Union{Real,Nothing} = nothing)

    _title = t!= nothing ? latexstring("⟨ϕ_0,ϕ( ",t," )⟩") : ti = L"⟨ϕ_0,ϕ(t)⟩"
    plt = plot(title = _title, size=(350,350),legend=false,framestyle=:origin)

    #Plot Unit circle
    θ = LinRange(0, 2π, 360)
    circle = (cos.(θ), sin.(θ))
    plot!(circle..., seriestype = [:shape], label="", fillcolor=:white, linewidth=2, fillalpha=0.0)

    # Plot arrow 
    plot!([0,real(Φ)], [0,imag(Φ)], arrow = true, color=:blue, linewidth=2, label="")

    plot!(draw_arrow=true)
    plot!(xlim=[-1.1,1.1],ylim=[-1.1,1.1],aspect_ratio=:equal)

    annotate!(-1.0,1.0,text(L"\mathbb{C}",20))

    plt
end

function plot_proba_l2Z(ψ::Vector, t::Union{Real,Nothing} = nothing; Nmax::Int, Nmin::Int = -Nmax)

    _title = t!= nothing ? latexstring("\text{Probabity of presence on each site of Z }",t) : _title = "Probabity of presence on each site at "
    plt = plot(title = _title, size=(800,350), legend = false)
    
    proba = abs.(ψ).^2
    M = 0.2 #max(proba[2:end]...) == 0.0 ? 1.0 : max(proba[2:end]...) 
    bar!(Nmin:Nmax, proba)

    plot!(ylim = [0.0,M])

    annotate!(0,M,text("Total probability presence = "*string(round(sum(proba),digits = 2)),10))
    annotate!(Nmin+2,M,text(L"\ell^2(\mathbb{Z})",20))

    plt
end

function plot_simu(ψ::Vector, t::Union{Real,Nothing} = nothing; info::NamedTuple)
    
    plt_C   = plot_compoC(ψ[1], t)
    plt_l2Z = plot_proba_l2Z(ψ[2:end], t; Nmax = info.Nmax, Nmin = info.Nmin)

    plot(plt_C,plt_l2Z, layout = (1,2), margin = 1Plots.cm, size=(1200,500))

    trace = ψ[1]
end