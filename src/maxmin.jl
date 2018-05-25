abstract type Downsampler{T, N} <: AbstractArray{T, N} end

struct MaxMin{T<:Number, S, A<:AbstractArray} <: Downsampler{Tuple{T, T}, 1}
    input::A
    binsize::Int
    len::Int
    input_len::Int
    function MaxMin{T,S,A}(input::A, binsize::Int) where
        {T<:Number, S<:Number, A<:AbstractVector{S}}
        input_len = length(input)
        len = cld(input_len, binsize)
        return new(input, binsize, len, input_len)
    end
    function MaxMin{T,S,A}(input::A, binsize::Int) where
        {T<:Number, S<:Number, A<:AbstractArray{S, 2}}
        @assert size(input, 1) == 2 "assumes max and min are on first dimension"
        input_len = size(input, 2)
        len = cld(input_len, binsize)
        return new(input, binsize, len, input_len)
    end
    function MaxMin{T,S,A}(input::A, binsize::Int) where
        {T<:Number, S<:NTuple{2, T}, A<:AbstractVector{S}}
        input_len = length(input)
        len = cld(input_len, binsize)
        return new(input, binsize, len, input_len)
    end
end
function MaxMin(a::A, n::Integer) where {T<:Number, A<:AbstractVector{T}}
    MaxMin{T, T, A}(a, convert(Int, n))
end
function MaxMin(a::A, n::Integer) where
    {T<:Number, S<:NTuple{2, T}, A<:AbstractVector{S}}
    MaxMin{T, S, A}(a, convert(Int, n))
end
function MaxMin(a::A, n::Integer) where
    {T<:Number, A<:AbstractArray{T, 2}}
    MaxMin{T, T, A}(a, convert(Int, n))
end

length(a::MaxMin) = a.len

size(a::MaxMin) = (length(a),)
binsize(a::MaxMin) = a.binsize

function getindex(a::MaxMin, i::Integer)
    (idx_start, idx_stop) = bin_bounds(i, a.binsize, a.input_len)
    return extrema_red(view(a.input, idx_start:idx_stop))
end
function getindex(a::M, i::Integer) where
    {T<: Number, A<:AbstractArray{T, 2}, M<:MaxMin{T, T, A}}
    (idx_start, idx_stop) = bin_bounds(i, a.binsize, a.input_len)
    return extrema_red(view(a.input, :, idx_start:idx_stop))
end

Base.IndexStyle(::Type{T}) where T<: MaxMin = IndexLinear()

setindex!(::MaxMin, ::Integer) = throw(ReadOnlyMemoryError())

function bin_bounds(a::MaxMin, ::Type{S} = Int) where {S}
    bnds = Vector{NTuple{2,S}}(length(a))
    bnds .= bin_bounds.(convert.(S, eachindex(a)), a.binsize, a.input_len)
    return bnds
end

function bin_bounds(i::Integer, a::MaxMin)
    bin_bounds(i, a.binsize, a.input_len)
end
