function prepare_cachefile(
    ::Type{<:Averager},
    input::AbstractArray{<:Any, N},
    dim::Integer = 1,
    reduce_dim::Bool = def_reduce(N)
) where {N}
    if reduce_dim
        ds = Averager(input, 10, dim)
        cachedims = (length(ds),)
    else
        ds = Averager(input, 10, dim)
        dims = collect(el_size(D, input, dim))
        push!(dims, length(ds))
        cachedims = (dims...)
    end
    return cachedims, ds
end
