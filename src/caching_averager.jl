function prepare_cachefile(
    ::Type{A},
    input::AbstractArray{<:Any, N},
    dim::Integer = N,
    reduce_dim::Bool = def_reduce(N)
) where {A<:Averager, N}
    ds = Averager(input, 10, dim)
    el_dims = el_size(A, input, dim, reduce_dim)

    ds_len = length(ds)
    if isempty(el_dims)
        cachedims = (ds_len,)
    else
        dim_vec = collect(el_dims)
        if dim_vec[end] == 1
            dim_vec[end] = ds_len
        else
            push!(dim_vec, ds_len)
        end
        cachedims = (dim_vec...)
    end
    return cachedims, ds
end
