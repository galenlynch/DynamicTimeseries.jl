struct Averager{
    E, T, N, A<:AbstractArray{T, N}
} <: Downsampler{E, 1}
    input::A
    dim::Int
    binsize::Int
    overlap::Int
    len::Int
    input_len::Int
    function Averager{E,T,N,A}(
        input::A, dim::Int, binsize::Int, overlap::Int = 0
    ) where {E,T,N,A<:AbstractArray{T,N}}
        dims = size(input)
        if dim <= 0 || dim > length(dims)
            throw(ArgumentError("dimension does not match input"))
        end
        input_len  = size(input, dim)
        if binsize < 0 || binsize > input_len
            throw(ArgumentError("binsize must be between 0 and input length"))
        end
        if overlap < 0 || overlap >= binsize
            throw(ArgumentError("overlap must be between 0 and binsize"))
        end
        # number of windows in signal
        # let N:= length of input
        # l := length of window
        # d := length of overlap
        # w := number of windows
        # Then, N = w * l - (w - 1) * d
        # So N = wl - wd + d = w(l-d) + d
        # w = (N - d) / (l - d)
        len = cld(input_len - overlap, binsize - overlap)
        return new(input, dim, binsize, overlap, len, input_len)
    end
end

function Averager(
    input::A, dim::Integer, binsize::Integer, overlap::Integer
) where {T,N,A<:AbstractArray{T,N}}
    S = div_type(T)
    E = Array{S, N}
    Averager{E,T,N,A}(
        input,
        convert(Int, dim),
        convert(Int, binsize),
        convert(Int, overlap)
    )
end

function Averager(
    input::A, binsize::Integer, overlap::Integer
) where {T,A<:AbstractVector{T}}
    S = div_type(T)
    Averager{S,T,1,A}(input, 1, convert(Int, binsize), convert(Int, overlap))
end

length(a::Averager) = a.len
size(a::Averager) = (a.len,)
Base.IndexStyle(::Type{T}) where T<:Averager = IndexLinear()
setindex!(::Averager, ::Integer) = throw(ReadOnlyMemoryError())

binsize(a::Averager) = a.binsize

function sliceview(a::Averager, i::Integer)
    (bb, be) = bin_bounds(i, a)
    idxes = make_slice_idx(a, bb, be)
    view(a.input, idxes...)
end

function getindex(a::Averager{E, <:Any, N, <:Any}, i::Integer) where {E,N}
    v = sliceview(a, i)
    mean(v, a.dim)::E
end

function getindex(a::Averager{E,<:Any,1,<:Any}, i::Integer) where E
    v = sliceview(a, i)
    mean(v)::E
end

function getindex(
    a::Averager{E,<:Any, N, <:Any}, idxes::OrdinalRange{<:Integer,<:Any}
) where {E,N}
    outdims = [size(a.input)...]
    outdims[a.dim] = length(idxes)
    out = E(outdims...)
    for i in idxes
        slice_idx = make_slice_idx(a, i)
        out[slice_idx...] = a[i]
    end
    out
end

function getindex(
    a::Averager{E,<:Any,1,<:Any}, idxes::OrdinalRange{<:Integer,<:Any}
) where {E}
    out = Vector{E}(length(idxes))
    for (i, idx) in enumerate(idxes)
        out[i] = a[idx]
    end
    out
end

function make_slice_idx(a::Averager, idx) 
    ndims = length(size(a.input))
    make_slice_idx(ndims, a.dim, idx)
end

function make_slice_idx(a::Averager, bb::T, be::T) where T<:Integer
    make_slice_idx(a, bb:be)
end

function bin_bounds(i::T, a::Averager) where {T<:Integer}
    start_offset = (i - one(T)) * T(a.binsize - a.overlap)
    end_idx = start_offset + T(a.binsize)
    return (start_offset + one(T), end_idx)
end
