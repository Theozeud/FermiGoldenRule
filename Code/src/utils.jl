# Utils function

apply(f::Base.Callable, v::AbstractVector) = f.(v)
apply(f::Base.Callable, v::Number) = f(v)
apply(f::Number, ::Nothing) = f

function vectorize(T::Tuple)
    vcat([e for e in T]...)
end

function vectorize2(T::Tuple)
    n = length(T)
    result = Vector{typeof(T[1][1])}(undef, n * length(T[1]))
    idx = 1
    for e in T
        for item in e
            result[idx] = item
            idx += 1
        end
    end
    return result
end
