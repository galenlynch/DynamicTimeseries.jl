function prepare_cachefile(
    ::Type{A},
    input::AbstractArray{<:Any, N},
    dim::Integer = 1,
    reduce_dim::Bool = def_reduce(N)
) where {A<:Averager, N}
    println("Preparing cachefile for averager")
    ds = Averager(input, 10, dim)
    el_dims = el_size(A, input, dim, reduce_dim)
    println("reduce_dim is ", reduce_dim)
    ds_len = length(ds)
    if isempty(el_dims)
        cachedims = (ds_len,)
    else
        dim_vec = collect(el_dims)
        push!(dim_vec, ds_len)
        cachedims = (dim_vec...)
    end
    return cachedims, ds
end
