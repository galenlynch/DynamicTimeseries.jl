struct MaxMin{T<:Number, W<:WindowedArray} <: Downsampler{Tuple{T, T}, 1}
    winput::W

    function MaxMin{T,W}(winput::W) where
        {T<:Number, W<:WindowedArray{<:Any,T,1,<:Any}}
        new(winput)
    end

    function MaxMin{T,W}(winput::W) where
        {T<:Number, S<:NTuple{2,T}, W<:WindowedArray{<:Any,S,<:Any,<:Any}}
        new(winput)
    end

    function MaxMin{T,W}(winput::W) where
        {T<:Number, W<:WindowedArray{<:Any,T,2,<:Any}}
        if size(winput.input, 1) != 2
            throw(ArgumentError("assumes max and min are on first dimension"))
        end
        new(winput)
    end
end

function MaxMin(a::W) where
    {T<:Number, S<:NTuple{2, T}, W<:WindowedArray{<:Any,S,1,<:Any}}
    MaxMin{T,W}(a)
end
function MaxMin(a::W) where
    {T<:Number, W<:WindowedArray{<:Any,T,<:Any,<:Any}}
    MaxMin{T,W}(a)
end

function MaxMin(input::AbstractArray{<:Number, 2}, binsize::Integer)
    winput = WindowedArray(input, binsize, 2)
    MaxMin(winput)
end
function MaxMin(input::AbstractVector, binsize::Integer)
    winput = WindowedArray(input, binsize)
    MaxMin(winput)
end

function concrete_type(
    ::Type{D}, ::Type{W}, args...
) where {D<:MaxMin, T, W<:WindowedArray{<:Any, T, <:Any, <:Any}}
    return MaxMin{T, W}
end

function eltype_preview(
    ::Type{D}, ::Type{<:AbstractVector{T}}
) where {D<:MaxMin, T<:Number}
    Tuple{T,T}
end

length(a::MaxMin) = length(a.winput)
size(a::MaxMin) = size(a.winput)

getindex(a::MaxMin, i::Integer) = extrema_red(a.winput[i])
setindex!(::MaxMin, ::Integer) = throw(ReadOnlyMemoryError())

Base.IndexStyle(::Type{M}) where {W,M<:MaxMin{<:Any,W}} = IndexStyle(W)

binsize(a::MaxMin) = binsize(a.winput)
bin_bounds(a::MaxMin, args...) = bin_bounds(a.winput, args...)
bin_bounds(i::Integer, a::MaxMin) = bin_bounds(i, a.winput)

downsamp_reduce(::Type{<:MaxMin}, ds::AbstractArray{<:Number, 2}) = extrema_red(ds)

function downsamp_reduce(
    ::Type{<:MaxMin}, ds::AbstractVector{<:NTuple{2, <:Number}}
)
    extrema_red(ds)
end

downsamp_reduce(::Type{<:MaxMin}, ds::AbstractVector{<:Number}) = extrema_red(ds)

function downsamp_reduce(
    ::Type{D}, ds::AbstractArray, ::AbstractVector, ::Integer = 0
) where D<:MaxMin
    return (downsamp_reduce(D, ds), 0)
end
