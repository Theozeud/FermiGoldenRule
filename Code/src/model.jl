# This files contains model and tool to  build model easily

function create_matrix(Hfun::Base.Callable, Infos_discr::Tuple)
    temp_arrays = base(Infos_discr)
    temp_arrays = [vectorize(Hfun(b...)) for b in temp_arrays]
    hcat(temp_arrays...)
end

