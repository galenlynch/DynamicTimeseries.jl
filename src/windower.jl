struct WindowedArray{E,T,N,A<:AbstractArray{T,N}}<: AbstractArray{E,1}
    input::A
    dim::Int
    binsize::Int
    overlap::Int
    len::Int
    input_len::Int
    function WindowedArray{E,T,N,A}(
        input::A, binsize::Int, dim::Int = 1, overlap::Int = 0
    ) where {E<:SubArray,T,N,A<:AbstractArray{T,N}}
        validate_windower_args(size(input), dim, binsize, overlap)
        input_len = size(input, dim)
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

function WindowedArray(
    input::A, binsize::Integer, dim::Integer = 1, overlap::Integer = 0
) where {T,N,A<:AbstractArray{T,N}}
    # make fake slice to determine SubArray type
    if dim < 1 || dim > N
        throw(ArgumentError("dimension out of bounds"))
    end
    si = make_slice_idx(N, dim, 1:0) # Empty slice to determine SubArray type
    v = view(input, si...)
    WindowedArray{typeof(v),T,N,A}(
        input,
        convert(Int, binsize),
        convert(Int, dim),
        convert(Int, overlap)
    )
end

function validate_windower_args(
    input_dims::NTuple{N,<:Integer}, dim::Integer, binsize::Integer, overlap::Integer
) where N
    if dim < 1 || dim > N
        throw(ArgumentError(string(
            "dim is ", dim,
            " but it must be between 1 and the dimension of the input, ", N
        )))
    end
    input_len  = input_dims[dim]
    if binsize < 1 || (input_len > 0 && binsize > input_len)
        throw(ArgumentError(string(
            "binsize is ", binsize,
            " but it must be between 1 and the input length, ", input_len
        )))
    end
    if overlap < 0 || overlap >= binsize
        throw(ArgumentError(string(
            "overlap is ", overlap,
            " but it must be between 0 and binsize, ", binsize
        )))
    end
end

length(a::WindowedArray) = a.len
size(a::WindowedArray) = (a.len,)
setindex!(::WindowedArray, ::Integer) = throw(ReadOnlyMemoryError())

Base.IndexStyle(::Type{T}) where T<:WindowedArray = IndexLinear()

binsize(a::WindowedArray) = a.binsize

function getindex(a::WindowedArray{E,<:Any,<:Any,<:Any}, i::Integer) where {E}
    sliceview(a, i)::E
end

function sliceview(a::WindowedArray, i::Integer)
    (bb, be) = bin_bounds(i, a)
    idxes = make_slice_idx(a, bb, be)
    view(a.input, idxes...)
end

function make_slice_idx(a::WindowedArray{<:Any,<:Any,N,<:Any}, idx)  where N
    make_slice_idx(N, a.dim, idx)
end

function make_slice_idx(a::WindowedArray, bb::T, be::T) where T<:Integer
    make_slice_idx(a, bb:be)
end

function bin_bounds(i::T, a::WindowedArray) where {T<:Integer}
    start_offset = (i - one(T)) * T(a.binsize - a.overlap)
    raw_end_idx = start_offset + T(a.binsize)
    end_idx = min(raw_end_idx, T(a.input_len))
    return (start_offset + one(T), end_idx)
end

function bin_bounds(a::WindowedArray, ::Type{T} = Int32) where T<:Integer
    bnds = Vector{NTuple{2,T}}(length(a))
    for i in eachindex(a)
        bnds[i] = bin_bounds(T(i), a)
    end
    bnds
end

basedata(a::WindowedArray) = a.input
