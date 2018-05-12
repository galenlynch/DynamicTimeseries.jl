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

function getindex(a::M, i::Integer) where
    {T<: Number, A<:AbstractVector, M<:MaxMin{T, T, A}}
    @boundscheck @assert checkbounds(Bool, a, i) "Index out of bounds"
    (idx_start, idx_stop) = bin_bounds(i, a.binsize, length(a.input))
    idx_range = idx_start:idx_stop
    return extrema(view(a.input, idx_range))
end
function getindex(a::M, i::Integer) where
    {T<: Number, S<:NTuple{2, T}, A<:AbstractVector, M<:MaxMin{T, S, A}}
    @boundscheck @assert checkbounds(Bool, a, i) "Index out of bounds"
    (idx_start, idx_stop) = bin_bounds(i, a.binsize, length(a.input))
    cmin = a.input[idx_start][1]
    cmax = a.input[idx_start][2]
    for i = (idx_start + 1):idx_stop
        cmin = min(cmin, a.input[i][1])
        cmax = max(cmax, a.input[i][2])
    end
    return (cmin, cmax)
end
function getindex(a::M, i::Integer) where
    {T<: Number, A<:AbstractArray{T, 2}, M<:MaxMin{T, T, A}}
    @boundscheck @assert checkbounds(Bool, a, i) "Index out of bounds"
    (idx_start, idx_stop) = bin_bounds(i, a.binsize, size(a.input, 2))
    idx_range = idx_start:idx_stop
    extremum = (
        minimum(view(a.input, 1, idx_range)),
        maximum(view(a.input, 2, idx_range))
    )
    return extremum
end

Base.IndexStyle(::Type{T}) where T<: MaxMin = IndexLinear()

setindex!(::MaxMin, ::Integer) = throw(ReadOnlyMemoryError())

bin_bounds(a::MaxMin) = bin_bounds.(eachindex(a), a.binsize, length(a.input))
bin_bounds(i::Integer, a::MaxMin) = bin_bounds(i, a.binsize, length(a.input))

"Must implement downsamp_req and duration"
abstract type DynamicDownsampler end

function downsamp_req(
    ds::DynamicDownsampler, x_start, x_end, reqpoints::AbstractFloat
)
    return downsamp_req(ds, x_start, x_end, floor(Int, reqpoints))
end

struct DynamicTs{A<:AbstractVector} <: DynamicDownsampler
    input::A
    fs::Float64
    offset::Float64
    function DynamicTs{A}(
        input::A,
        fs::Float64,
        offset::Float64
    ) where {A<:AbstractArray}
        return new(input, fs, offset)
    end
end
function DynamicTs(
    input::A,
    fs::Real,
    offset::Real = 0
) where {A <: AbstractVector}
    return DynamicTs{A}(input, convert(Float64, fs), convert(Float64, offset))
end

function downsamp_req(dts::DynamicTs, x_start, x_end, reqpoints::Integer)
    nin = length(dts.input)

    # Find bounding indices
    i_begin = clip_ndx(t_to_ndx(x_start, dts.fs, dts.offset), nin)
    i_end = clip_ndx(t_to_ndx(x_end, dts.fs, dts.offset), nin)

    # Calculate binsize
    npt = n_ndx(i_begin, i_end)
    binsize = cld(npt, reqpoints)

    was_downsampled = reqpoints > npt

    # Clip data
    subview = view(dts.input, i_begin:i_end)

    # Create downsampled data
    ys = MaxMin(subview, binsize)
    bin_start_x = ndx_to_t(i_begin, dts.fs, dts.offset)
    xs = ndx_to_t(bin_center(bin_bounds(ys)), dts.fs, bin_start_x)

    return (xs, ys, was_downsampled)
end

duration(d::DynamicTs) = duration(length(d.input), d.fs, d.offset)

struct CachingDynamicTs{S<:Number, A<:AbstractVector{S}} <: DynamicDownsampler
    input::A
    fs::Float64
    offset::Float64
    cachearrays::Vector{Array{S, 2}}
    function CachingDynamicTs{S,A}(
        input::A,
        fs::Float64,
        offset::Float64,
        cachearrays::Vector{Array{S, 2}},
    ) where {S<:Number, A<:AbstractVector{S}}
        return new(input, fs, offset, cachearrays)
    end
end
function CachingDynamicTs(
    input::A,
    fs::Real,
    offset::Real,
    cachearrs::Vector{Array{S, 2}},
) where {S<:Number, A<:AbstractVector{S}}
    return CachingDynamicTs{S,A}(
        input,
        convert(Float64, fs),
        convert(Float64, offset),
        cachearrs
    )
end
function CachingDynamicTs(
    input::A,
    fs::Real,
    offset::Real = 0,
    sizehint::Integer = 300, # x dimension in pixels of a small window?
    autoclean::Bool = true
) where {S<:Number, A<:AbstractVector{S}}
    (cachepaths, cachelengths) = write_cache_files(input, sizehint, autoclean)
    cachearrs = open_cache_files(S, cachelengths, cachepaths, false)
    CachingDynamicTs(input, fs, offset, cachearrs)
end

function downsamp_req(dts::CachingDynamicTs, x_start, x_end, reqpts::Integer)
    nin = length(dts.input)
    i_begin = clip_ndx(t_to_ndx(x_start, dts.fs, dts.offset), nin)
    i_end = clip_ndx(t_to_ndx(x_end, dts.fs, dts.offset), nin)
    nbase = n_ndx(i_begin, i_end)
    ncache = length(dts.cachearrays)
    decno = min(floor(Int, log10(nbase / reqpts)), ncache)
    if decno >= 1
        was_downsampled = true
        di_begin = decade_ndx_conversion(i_begin, decno)
        di_end = decade_ndx_conversion(i_end, decno)
        subview = view(dts.cachearrays[decno], 1:2, di_begin:di_end)

        nsubsamp = n_ndx(di_begin, di_end)
        binsize = cld(nsubsamp, reqpts)
        ys = MaxMin(subview, binsize)

        dec_idxes = bin_center(bin_bounds(ys)) + di_begin - 1
        base_idxes = dec_ndx_to_ndx(dec_idxes, decno)
        xs = ndx_to_t(base_idxes, dts.fs, dts.offset)
    else
        binsize = ceil(Int, nbase / reqpts)
        was_downsampled = binsize > 1
        i_range = i_begin:i_end

        # Clip data
        subview = view(dts.input, i_range)
        ys = MaxMin(subview, binsize)

        bin_start_x = ndx_to_t(i_begin, dts.fs, dts.offset)
        xs = ndx_to_t(bin_center(bin_bounds(ys)), dts.fs, bin_start_x)
    end
    return (xs, ys, was_downsampled)
end

duration(dts::CachingDynamicTs) = duration(length(dts.input), dts.fs, dts.offset)

decade_ndx_conversion(i::Integer, dec::Integer) = fld(i - 1, 10 ^ dec) + 1
ndx_to_dec_ndx(args...) = decade_ndx_conversion(args...)
dec_ndx_to_ndx(i::Real, dec::Integer) = bin_center(i, 10 ^ dec)
dec_ndx_to_ndx(is::AbstractVector, dec::Integer) = bin_center.(is, 10 ^ dec)

struct MappedDynamicDownsampler{D<:DynamicDownsampler} <: DynamicDownsampler
    downsampler::D
    fmap::Function
end

function downsamp_req(mdds::D, xb, xe, npts::Integer) where
    {D<:MappedDynamicDownsampler}
    (xs, ys, was_downsampled) = downsamp_req(mdds.downsampler, xb, xe, npts)
    mys = mdds.fmap(ys)
    return (xs, mys, was_downsampled)
end

duration(d::MappedDynamicDownsampler) = duration(d.downsampler)

function shift_extrema(shift, ys::T) where {S, T<:NTuple{2,S}}
    (ys[1] + shift, ys[2] + shift)
end
function shift_extrema(shift, ys::A) where
    {T, S<:NTuple{2,T}, A<:AbstractVector{S}}
    shifted = similar(ys)
    shift_extrema!(shift, shifted, ys)
    return shifted
end

function shift_extrema!(
    shift,
    dest::B,
    ys::A
) where {T, S<:NTuple{2,T}, A<:AbstractVector{S}, B<:AbstractVector{S}}
    for (i, tup) in enumerate(ys)
        dest[i] = shift_extrema(shift, ys[i])
    end
end

make_shifter(shift) = (x) -> shift_extrema(shift, x)
