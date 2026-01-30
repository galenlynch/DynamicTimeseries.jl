# D (bool) determines whether D is reduced
struct Averager{E,W<:WindowedArray,D} <: Downsampler{E,1}
    winput::W
    function Averager{E,W,D}(winput::W) where {E,T,N,W<:WindowedArray{<:Any,T,N,<:Any},D}
        isa(D, Bool) || throw(ArgumentError{"parameter D must be true or false"})
        exp_type = eltype_preview(Averager, T, N, D)
        if E != exp_type
            throw(ArgumentError(string("Expected element type ", E, " but got ", exp_type)))
        end
        new(winput)
    end
end

function Averager(
    winput::W,
    reduce_dim::Bool = def_reduce(N),
) where {T,N,W<:WindowedArray{<:Any,T,N,<:Any}}
    E = eltype_preview(Averager, T, N, reduce_dim)
    Averager{E,W,reduce_dim}(winput)
end

def_reduce(N::Integer) = N > 1 ? false : true

function Averager(input::AbstractArray, binsize::Integer, overlap::Integer = 0, args...)
    winput = WindowedArray(input, binsize, overlap)
    Averager(winput, args...)
end

length(a::Averager) = length(a.winput)
size(a::Averager) = size(a.winput)

function concrete_type(
    ::Type{D},
    ::Type{W},
    reduce_dim::Bool = def_reduce(N),
) where {D<:Averager,T,N,W<:WindowedArray{<:Any,T,N,<:Any}}
    E = eltype_preview(D, T, N, reduce_dim)
    return Averager{E,W,reduce_dim}
end

function eltype_preview(
    ::Type{<:Averager},
    ::Type{T},
    N::Integer,
    reduce_dim::Bool = def_reduce(N),
) where {T}
    S = div_type(T)
    reduce_dim ? S : Array{S,N - 1}
end

function eltype_preview(::Type{D}, ::Type{<:AbstractArray{T,N}}) where {D<:Averager,T,N}
    eltype_preview(D, T, N)
end

function arr_eltype_preview(::Type{<:Averager}, ::Type{A}) where {A<:AbstractArray}
    T = arr_bit_eltype(A)
    div_type(T)
end

function el_size(
    ::Type{<:Averager},
    a::AbstractArray{E,N},
    dim::Integer,
    dim_red::Bool = def_reduce(N),
) where {E,N}
    dim_red ? reduce_dim_size(a) : noreduce_dim_size(a, dim)
end

function el_size(
    ::Type{<:Averager{<:Any,<:Any,true}},
    a::AbstractArray,
    ::Integer,
    ::Bool = false,
)
    noreduce_dim_size(a, i)
end

function el_size(
    ::Type{<:Averager{<:Any,<:Any,false}},
    a::AbstractArray,
    dim::Integer,
    ::Bool = false,
)
    reduce_dim_size(a)
end

function noreduce_dim_size(a::AbstractArray, dim::Integer)
    el_dims = reduce_dim_size(a)
    dims = collect(size(a))
    dims[dim] = 1
    if ! isempty(el_dims)
        el_dim_vec = collect(el_dims)
        append!(el_dim_vec, dims)
        dims = el_dim_vec
    end
    (dims...,)
end

function reduce_dim_size(a::AbstractArray{<:AbstractArray,<:Any})
    # assumes elements are all the same size
    if isempty(a)
        throw(ArgumentError("Empty array"))
    end
    nested_arr_size(a[1])
end
reduce_dim_size(a::AbstractArray{<:Number,<:Any}) = ()
reduce_dim_size(a::AbstractArray, ::Integer) = reduce_dim_size(a)

function nested_arr_size(a::AbstractArray{<:AbstractArray,<:Any})
    # assumes elements are all the same size
    if ! isempty(a)
        el_dims = collect(nested_arr_size(a[1]))
        append!(el_dims, collect(size(a)))
    else
        throw(ArgumentError("Empty array"))
    end
    return (el_dims...,)
end
nested_arr_size(a::AbstractArray{<:Number,<:Any}) = size(a)

function getindex(a::Averager{E,<:Any,true}, i::Integer) where {E}
    v = a.winput[i]
    mean(v)::E
end

"Method meant for when N > 1 in WindowedArray"
function getindex(
    a::Averager{E,<:WindowedArray{<:Any,<:Any,N,<:Any},false},
    i::Integer,
) where {E,N}
    v = a.winput[i]
    dropdims(mean(v, dims = N); dims = N)::E
end
"Method meant for when N > 1 in WindowedArray"
function getindex(
    a::Averager{E,<:WindowedArray{<:Any,<:Any,1,<:Any},false},
    i::Integer,
) where {E}
    mean(a.winput[i])::E
end

setindex!(::Averager, ::Integer) = throw(ReadOnlyMemoryError())

Base.IndexStyle(::Type{T}) where {W,T<:Averager{<:Any,W}} = IndexStyle(W)
binsize(a::Averager) = binsize(a.winput)

bin_bounds(a::Averager, args...) = bin_bounds(a.winput, args...)
bin_bounds(i::Integer, a::Averager) = bin_bounds(i, a.winput)

function downsamp_reduce(
    ::Type{<:Averager{<:Any,<:Any,false}},
    ds::AbstractArray{<:Any,N},
    weigths::AbstractVector{<:Number},
) where {N}
    if size(ds, N) != length(weights)
        throw(
            ArgumentError(
                string(
                    "ds (size $(size(ds))) ",
                    "and weights (size $(size(weights))) ",
                    "are not the same size on dim $dim",
                ),
            ),
        )
    end
    total_weight = sum(weights)
    reduced = weighted_mean_dim(ds, weights, N, total_weight)
    return (reduced, total_weight)
end

function downsamp_reduce(
    ::Type{<:Averager{<:Any,<:Any,true}},
    ds::AbstractArray{<:Any,N},
    weights::AbstractArray{<:Number,N},
    ::Integer = 0,
) where {N}
    if size(ds) != size(weights)
        throw(
            ArgumentError(
                string(
                    "ds (size $(size(ds))) ",
                    "and weights (size $(size(weights))) ",
                    "are not the same size",
                ),
            ),
        )
    end
    total_weight = sum(weights)
    reduced = weighted_mean(ds, weights, total_weight)
    return (reduced, total_weight)
end

function downsamp_reduce_cache(
    ::Type{<:Averager{<:Any,<:Any,true}},
    ds::AbstractArray{<:Any,N},
    weights::AbstractVector{<:Number},
) where {N}
    if size(ds, N) != length(weights)
        throw(ArgumentError("last dim of ds does not match weights"))
    end
    total_weight = sum(weights)
    reduced = weighted_mean_dim(ds, weights, N, total_weight)
    squeezed = dropdims(reduced, dims = N)
    return (squeezed, total_weight)
end

function downsamp_reduce_cache(
    ::Type{A},
    ds::AbstractVector,
    weights::AbstractVector{<:Number},
) where {A<:Averager{<:Any,<:Any,true}}
    downsamp_reduce(A, ds, weights)
end
