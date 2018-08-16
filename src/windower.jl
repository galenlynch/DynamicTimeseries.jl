struct WindowedArray{E,T,N,A<:AbstractArray{T,N}}<: AbstractArray{E,1}
    input::A
    binsize::Int
    overlap::Int
    len::Int
    input_len::Int
    function WindowedArray{E,T,N,A}(
        input::A, binsize::Int, overlap::Int = 0
    ) where {E<:SubArray,T,N,A<:AbstractArray{T,N}}
        validate_windower_args(size(input), binsize, overlap)
        input_len = size(input, N)
        # number of windows in signal
        # let N:= length of input
        # l := length of window
        # d := length of overlap
        # w := number of windows
        # Then, N = w * l - (w - 1) * d
        # So N = wl - wd + d = w(l-d) + d
        # w = (N - d) / (l - d)
        len = cld(input_len - overlap, binsize - overlap)
        return new(input, binsize, overlap, len, input_len)
    end
end

function WindowedArray(
    input::A, binsize::Integer, overlap::Integer = 0
) where {T,N,A<:AbstractArray{T,N}}
    # make fake slice to determine SubArray type
    v = view_trailing_slice(input, 1:0)
    WindowedArray{typeof(v),T,N,A}(
        input,
        convert(Int, binsize),
        convert(Int, overlap)
    )
end

function validate_windower_args(
    input_dims::NTuple{N,<:Integer}, binsize::Integer, overlap::Integer
) where N
    input_len  = input_dims[N]
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
    view_trailing_slice(a.input, bb:be)
end

function bin_bounds(i::T, a::WindowedArray) where {T<:Integer}
    start_offset = (i - one(T)) * T(a.binsize - a.overlap)
    raw_end_idx = start_offset + T(a.binsize)
    end_idx = min(raw_end_idx, T(a.input_len))
    return (start_offset + one(T), end_idx)
end

function bin_bounds(a::WindowedArray, ::Type{T} = Int32) where T<:Integer
    @compat bnds = Vector{NTuple{2,T}}(undef, length(a))
    for i in eachindex(a)
        bnds[i] = bin_bounds(T(i), a)
    end
    bnds
end

basedata(a::WindowedArray) = a.input
