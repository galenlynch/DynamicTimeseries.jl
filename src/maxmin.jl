abstract type Downsampler{T, N} <: AbstractArray{T, N} end

struct MaxMin{T<:Number, S, A<:AbstractArray} <: Downsampler{Tuple{T, T}, 1}
    input::A
    binsize::Int
    function MaxMin{T,S,A}(input::A, binsize::Integer) where
        {T<:Number, S<:Number, A<:AbstractVector{S}}
        return new(input, convert(Int, binsize))
    end
    function MaxMin{T,S,A}(input::A, binsize::Integer) where
        {T<:Number, S<:Number, A<:AbstractArray{S, 2}}
        @assert size(input, 1) == 2 "assumes max and min are on first dimension"
        return new(input, convert(Int, binsize))
    end
    function MaxMin{T,S,A}(input::A, binsize::Integer) where
        {T<:Number, S<:NTuple{2, T}, A<:AbstractVector{S}}
        return new(input, convert(Int, binsize))
    end
end
function MaxMin(a::A, n::Integer) where {T<:Number, A<:AbstractVector{T}}
    MaxMin{T, T, A}(a, n)
end
function MaxMin(a::A, n::Integer) where
    {T<:Number, S<:NTuple{2, T}, A<:AbstractVector{S}}
    MaxMin{T, S, A}(a, n)
end
function MaxMin(a::A, n::Integer) where
    {T<:Number, A<:AbstractArray{T, 2}}
    MaxMin{T, T, A}(a, n)
end

function length(a::M) where {T, S, A<:AbstractVector, M<:MaxMin{T, S, A}}
    cld(length(a.input), a.binsize)
end
function length(a::M) where {T, S, A<:AbstractArray{T, 2}, M<:MaxMin{T, S, A}}
    cld(size(a.input, 2), a.binsize)
end

size(a::MaxMin) = (length(a),)
binsize(a::MaxMin) = a.binsize

function getindex(a::MaxMin, i::Integer)
    (idx_start, idx_stop) = bin_bounds(i, a.binsize, length(a.input))
    return extrema_red(view(a.input, idx_start:idx_stop))
end
function getindex(a::M, i::Integer) where
    {T<: Number, A<:AbstractArray{T, 2}, M<:MaxMin{T, T, A}}
    (idx_start, idx_stop) = bin_bounds(i, a.binsize, size(a.input, 2))
    return extrema_red(view(a.input, :, idx_start:idx_stop))
end

Base.IndexStyle(::Type{T}) where T<: MaxMin = IndexLinear()

setindex!(::MaxMin, ::Integer) = throw(ReadOnlyMemoryError())

function bin_bounds(a::MaxMin, ::Type{S} = Int) where {S}
    bnds = Vector{NTuple{2,S}}(length(a))
    bnds .= bin_bounds.(convert.(S, eachindex(a)), a.binsize, length(a.input))
    return bnds
end

function bin_bounds(i::Integer, a::MaxMin, ::Type{T}) where {T}
    bin_bounds(i, a.binsize, length(a.input))
end
