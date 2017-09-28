abstract type Downsampler{T, N} <: AbstractArray{T, N} end

immutable MaxMin{T<:Number, A<:AbstractVector{T}} <: Downsampler{Tuple{T, T}, 1}
    input::A
    binsize::Int
    function MaxMin{T, A}(input::A, binsize::Integer) where {T<:Number, A<:AbstractVector{T}}
        return new(input, convert(Int, binsize))
    end
end
MaxMin(a::A, n::Integer) where {T<:Number, A<:AbstractVector{T}} = MaxMin{T, A}(a, n)

length(a::MaxMin) = cld(length(a.input), a.binsize)

size(a::MaxMin) = (length(a),)
binsize(a::MaxMin) = a.binsize

function getindex(a::MaxMin, i::Int)
    @boundscheck @assert checkbounds(Bool, a, i) "Index out of bounds"
    (idx_start, idx_stop) = bin_bounds(i, a.binsize, length(a.input))
    idx_range = idx_start:idx_stop
    return extrema(view(a.input, idx_range))
end

Base.IndexStyle(::Type{T}) where T<: MaxMin = IndexLinear()

setindex!(::MaxMin, ::Int) = throw(ReadOnlyMemoryError())

bin_bounds(a::MaxMin) = bin_bounds.(eachindex(a), a.binsize, length(a.input))
bin_bounds(i::Integer, a::MaxMin) = bin_bounds(i, a.binsize, length(a.input))

immutable DynamicTs{T<:Number, A<:AbstractVector{T}}
    input::A
    fs::Float64
    offset::Float64
    function DynamicTs{T, A}(
        input::A,
        fs::Float64,
        offset::Float64
    ) where {T<:Number, A<:AbstractVector{T}}
        return new(input, fs, offset)
    end
end
function DynamicTs(
    input::A,
    fs::Real,
    offset::Real
) where {T <: Number, A <: AbstractVector{T}}
    return DynamicTs{T, A}(input, convert(Float64, fs), convert(Float64, offset))
end
DynamicTs(input::OEContArray, offset::Real = 0) = DynamicTs(input, input.contfile.header.samplerate, offset)

function downsamp_req(dts::DynamicTs, x_start::Real, x_end::Real, maxpoints::Integer)
    # Find bounding indices
    i_start = x_to_ndx(x_start, dts.fs, dts.offset)
    i_end = x_to_ndx(x_end, dts.fs, dts.offset)

    # Calculate binsize
    npt = n_ndx(i_start, i_end)
    binsize = convert(Int, ceil(npt / maxpoints))

    # Clip data
    subview = view(dts.input, i_start:i_end)

    # Create downsampled data
    ys = MaxMin(subview, binsize)
    bin_start_x = ndx_to_x(i_start, dts.fs, dts.offset)
    xs = ndx_to_x(bin_center(bin_bounds(ys)), dts.fs, bin_start_x)

    return (xs, ys)
end
function downsamp_req(dts::DynamicTs, x_start::Real, x_end::Real, maxpoints::AbstractFloat)
    return downsamp_req(dts, x_start, x_end, convert(Int, floor(maxpoints)))
end
