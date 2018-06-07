function prepare_cachefile(
    ::Type{D}, input::A, dim::Integer = 1
) where {D<:Averager, A<:AbstractVector}
    S = eltype_preview(D, A)
    ds = Averager(input, 10, dim)
    return (length(ds),), S, ds
end

function prepare_cachefile(
    ::Type{D}, input::A, dim::Integer = 1
) where {D<:Averager, A<:AbstractArray}
    S = eltype_preview(D, A)
    N = array_n(S)
    ds = Averager(input, 10, dim)
    dims = collect(el_size(D, input, dim))
    dims[dim] = 1
    return ((dims...), S, ds)
end

array_n(::Type{<:AbstractArray{<:Any, N}}) where {N} = N
