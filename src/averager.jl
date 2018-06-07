# D (bool) determines whether D is reduced
struct Averager{E,W<:WindowedArray,D} <: Downsampler{E,1}
    winput::W
    function Averager{E,W,D}(
        winput::W
    ) where {E,T,N,W<:WindowedArray{<:Any,T,N,<:Any}, D}
        isa(D, Bool) || throw(ArgumentError{"parameter D must be true or false"})
        S = div_type(T)
        exp_type = D ? S : Array{S,N}
        if E != exp_type
            throw(ArgumentError(
                "Expected element type ", E, " but got ", exp_type
            ))
        end
        new(winput)
    end
end

function Averager(
    winput::W, reduce_dim::Bool = def_reduce(N)
) where {T,N,W<:WindowedArray{<:Any,T,N,<:Any}}
    S = div_type(T)
    E = reduce_dim ? S : Array{S,N}
    Averager{E,W,reduce_dim}(winput)
end

def_reduce(N::Integer) = N > 1 ? false : true

function Averager(
    input::AbstractArray,
    binsize::Integer,
    dim::Integer = 1,
    overlap::Integer=0,
    args...
)
    winput = WindowedArray(input, binsize, dim, overlap)
    Averager(winput, args...)
end

length(a::Averager) = length(a.winput)
size(a::Averager) = size(a.winput)

function concrete_type(
    ::Type{D}, ::Type{W}, reduce_dim::Bool = def_reduce(N)
) where {D<:Averager, T, N, W<:WindowedArray{<:Any, T, N, <:Any}}
    E = eltype_preview(D, T, N, reduce_dim)
    return Averager{E, W, reduce_dim}
end

function eltype_preview(
    ::Type{<:Averager}, ::Type{T}, N::Integer, reduce_dim::Bool = def_reduce(N)
) where {T<:Number}
    S = div_type(T)
    reduce_dim ? S : Array{S,N}
end

function eltype_preview(
    ::Type{D},
    ::Type{<:AbstractArray{T,N}}
) where {D<:Averager, T<:Number, N}
    eltype_preview(D, T, N)
end

el_size(::Type{<:Averager{<:Any,<:Any,true}}, a::AbstractArray, ::Integer) = ()

function el_size(::Type{<:Averager{<:Any,<:Any,false}}, a::AbstractArray, ::Integer)
    dims = collect(size(a))
    dims[dim] = 1
    (dims...)
end

function getindex(
    a::Averager{E, <:Any, true}, i::Integer
) where E
    v = a.winput[i]
    mean(v)::E
end

"Method meant for when N > 1 in WindowedArray"
function getindex(
    a::Averager{E, <:WindowedArray, false}, i::Integer
) where E
    v = a.winput[i]
    mean(v, a.winput.dim)::E
end
setindex!(::Averager, ::Integer) = throw(ReadOnlyMemoryError())

Base.IndexStyle(::Type{T}) where {W,T<:Averager{<:Any,W}} = IndexStyle(W)
binsize(a::Averager) = binsize(a.winput)

bin_bounds(a::Averager, args...) = bin_bounds(a.winput, args...)
bin_bounds(i::Integer, a::Averager) = bin_bounds(i, a.winput)

function downsamp_reduce(
    ::Type{<:Averager{<:Any, <:Any, false}},
    ds::AbstractArray{E, N},
    weigths::AbstractVector{T},
    dim::Integer = N
) where {E<:Number, N, T<:Number}
    if size(ds, dim) != length(weights)
        throw(ArgumentError("ds and weights are not the same size"))
    end
    total_weight = sum(weights)
    reduced = weighted_mean_dim(ds, weights, dim, total_weight)
    return (reduced, total_weight)
end

function downsamp_reduce(
    ::Type{<:Averager{<:Any, <:Any, true}},
    ds::AbstractVector{<:Any},
    weights::AbstractVector{<:Number},
    ::Integer = 0
)
    if length(ds) != length(weights)
        throw(ArgumentError("ds and weights are not the same size"))
    end
    total_weight = sum(weights)
    reduced = weighted_mean(ds, weights, total_weight)
    return (reduced, total_weight)
end
