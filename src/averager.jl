struct Averager{E,W<:WindowedArray} <: Downsampler{E,1}
    winput::W
end

function Averager(winput::W) where {T,N,W<:WindowedArray{<:Any,T,N,<:Any}}
    S = div_type(T)
    E = Array{S,N}
    Averager{E,W}(winput)
end
function Averager(winput::W) where {T,W<:WindowedArray{<:Any,T,1,<:Any}}
    S = div_type(T)
    Averager{S,W}(winput)
end

function Averager(
    input::AbstractArray, binsize::Integer, dim::Integer = 1, overlap::Integer=0
)
    winput = WindowedArray(input, binsize, dim, overlap)
    Averager(winput)
end

length(a::Averager) = length(a.winput)
size(a::Averager) = size(a.winput)


function getindex(
    a::Averager{E,<:WindowedArray{<:Any,<:Any,1,<:Any}}, i::Integer
) where E
    v = a.winput[i]
    mean(v)::E
end

"Method meant for when N > 1 in WindowedArray"
function getindex(a::Averager{E, <:WindowedArray}, i::Integer) where E
    v = a.winput[i]
    mean(v, a.winput.dim)::E
end
setindex!(::Averager, ::Integer) = throw(ReadOnlyMemoryError())

Base.IndexStyle(::Type{T}) where {W,T<:Averager{<:Any,W}} = IndexStyle(W)
binsize(a::Averager) = binsize(a.winput)

bin_bounds(a::Averager, args...) = bin_bounds(a.winput, args...)
bin_bounds(i::Integer, a::Averager) = bin_bounds(i, a.winput)

