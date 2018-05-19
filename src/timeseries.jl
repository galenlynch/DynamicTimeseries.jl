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

bin_bounds(a::MaxMin) = bin_bounds.(eachindex(a), a.binsize, length(a.input))
bin_bounds(i::Integer, a::MaxMin) = bin_bounds(i, a.binsize, length(a.input))

"Must implement downsamp_req and duration and extrema"
abstract type DynamicDownsampler{E} end

function extrema(d::DynamicDownsampler)
    (_, ys, _) = downsamp_req(d, duration(d)..., 1)
    extrema_red(ys)
end

function downsamp_req(
    ds::DynamicDownsampler, x_start, x_end, reqpoints::AbstractFloat
)
    return downsamp_req(ds, x_start, x_end, floor(Int, reqpoints))
end

struct DynamicTs{S<:Number, A<:AbstractVector{S}} <: DynamicDownsampler{NTuple{2, S}}
    input::A
    fs::Float64
    offset::Float64
    function DynamicTs{S, A}(
        input::A,
        fs::Float64,
        offset::Float64
    ) where {S<:Number, A<:AbstractVector{S}}
        return new(input, fs, offset)
    end
end
function DynamicTs(
    input::A,
    fs::Real,
    offset::Real = 0
) where {S<:Number, A <: AbstractVector{S}}
    return DynamicTs{S,A}(input, convert(Float64, fs), convert(Float64, offset))
end

function downsamp_range_check(dts::DynamicDownsampler, x_start, x_end)
    x_start <= x_end || throw(ArgumentError("x_start must be before x_end"))
    nin = length(dts.input)
    # Find bounding indices
    i_begin = t_to_ndx(x_start, dts.fs, dts.offset)
    i_end = t_to_ndx(x_end, dts.fs, dts.offset)
    in_range = i_begin <= nin && i_end >= 1
    i_begin_clipped = clip_ndx(i_begin, nin)
    i_end_clipped = clip_ndx(i_end, nin)
    return (in_range, i_begin_clipped, i_end_clipped)
end

function downsamp_req(
    dts::DynamicTs{S,A}, x_start, x_end, reqpoints::Integer
) where {S, A}
    (in_range, i_begin, i_end) = downsamp_range_check(dts, x_start, x_end)
    if ! in_range
        return (Vector{Float64}(), Vector{NTuple{2,S}}(), false)
    end
    # Calculate binsize
    npt = n_ndx(i_begin, i_end)
    binsize = max(fld(npt, reqpoints), 1)

    was_downsampled = binsize > 1

    # Clip data
    subview = view(dts.input, i_begin:i_end)

    # Create downsampled data
    mm = MaxMin(subview, binsize)
    bin_start_x = ndx_to_t(i_begin, dts.fs, dts.offset)
    xs = ndx_to_t(bin_center(bin_bounds(mm)), dts.fs, bin_start_x)
    ys = collect(mm)
    return (xs::Vector{Float64}, ys::Vector{NTuple{2,S}}, was_downsampled)
end

duration(d::DynamicTs) = duration(length(d.input), d.fs, d.offset)

struct CachingDynamicTs{S<:Number, A<:AbstractVector{S}} <: DynamicDownsampler{NTuple{2, S}}
    input::A
    fs::Float64
    offset::Float64
    cachearrays::Vector{Array{S, 2}}
    cachepaths::Vector{String}
    function CachingDynamicTs{S,A}(
        input::A,
        fs::Float64,
        offset::Float64,
        cachearrays::Vector{Array{S, 2}},
        cachepaths::Vector{String}
    ) where {S<:Number, A<:AbstractVector{S}}
        if length(cachepaths) != length(cachearrays)
            throw(ArgumentError("cachearrays must be same length as cachepaths"))
        end
        return new(input, fs, offset, cachearrays, cachepaths)
    end
end

function CachingDynamicTs(
    input::A,
    fs::Real,
    offset::Real,
    cachearrs::Vector{Array{S, 2}},
    cachepaths::Vector{String}
) where {S<:Number, A<:AbstractVector{S}}
    return CachingDynamicTs{S,A}(
        input,
        convert(Float64, fs),
        convert(Float64, offset),
        cachearrs,
        cachepaths
    )
end

function CachingDynamicTs(
    input::A,
    fs::Real,
    offset::Real,
    cachepaths::AbstractArray{<:AbstractString},
    cachelengths::AbstractArray{<:Integer},
    autoclean::Bool = false;
    checkfiles::Bool = true
) where {S<:Number, A<:AbstractVector{S}}
    if checkfiles
        p = sortperm(cachelengths; rev=true)
        cachelengths = cachelengths[p]
        cachepaths = cachepaths[p]
        last_len = length(input)
        for l in cachelengths
            if cld(last_len, l) != 10
                error(
"last cache was $last_len, but this cache is $l which is not a factor of 10"
                )
            end
            last_len = l
        end
    end
    cachearrs = open_cache_files(S, cachepaths, cachelengths, autoclean)
    CachingDynamicTs(input, fs, offset, cachearrs, cachepaths)
end

function CachingDynamicTs(
    input::AbstractVector{<:Number},
    fs::Real,
    offset::Real = 0,
    sizehint::Integer = 70, # x dimension in pixels of a small window?
    autoclean::Bool = true
)
    (cachepaths, cachelengths) = write_cache_files(input, sizehint, autoclean)
    CachingDynamicTs(
        input, fs, offset, cachepaths, cachelengths, false; checkfiles=false
    )
end

dec_ndx_greater(i, dec) = cld(i - 1, 10 ^ dec) + 1
dec_ndx_lesser(i, dec) = fld(i, 10 ^ dec)
dec_ndx_2_base_start(i, dec) = (i - 1) * 10 ^ dec + 1
dec_ndx_2_base_end(i, dec) = i * 10 ^ dec

"Recursively find extrema using cached arrays"
function extrema_range(dts::CachingDynamicTs{S, A}, i_start::Integer, i_stop::Integer) where {S, A}
    nbase = n_ndx(i_start, i_stop)
    decno = min(floor(Int, log10(nbase)), length(dts.cachearrays))
    if decno > 0
        cached_extrema = Vector{NTuple{2, S}}(3) # Store
        n_input = length(dts.input)

        # Find which cached indexes are completely contained in this slice
        contained_first = dec_ndx_greater(i_start, decno)
        if n_input == i_stop
            contained_last = size(dts.cachearrays[decno], 2)
        else
            contained_last = dec_ndx_lesser(i_stop, decno)
        end

        # Reduce contained cache values
        contained_view = view(dts.cachearrays[decno], :, contained_first:contained_last)
        if size(contained_view, 2) > 0
            extrema_i = 2
            cached_extrema[1] = extrema_red(contained_view)
            remainder_left_stop = dec_ndx_2_base_start(contained_first, decno) - 1
            remainder_right_start = dec_ndx_2_base_end(contained_last, decno) + 1
        else
            extrema_i = 1
            pop!(cached_extrema)
            # Split indicies in a cache-friendly manner
            next_dec = decno - 1
            remainder_left_stop = (10 ^ next_dec) * decade_ndx_conversion(
                i_start + fld(nbase, 2), next_dec) # Should be at the end of a cache
            remainder_right_start = remainder_left_stop + 1
        end

        # Recurse on 'left' remainder
        if i_start <= remainder_left_stop
            cached_extrema[extrema_i] = extrema_range(dts, i_start, remainder_left_stop)
            extrema_i += 1
        else
            pop!(cached_extrema)
        end

        # Recurse on 'right' remainder
        if i_stop >= remainder_right_start
            cached_extrema[extrema_i] = extrema_range(dts, remainder_right_start, i_stop)
            extrema_i += 1
        else
            pop!(cached_extrema)
        end

        out = extrema_red(cached_extrema)
    else # We're at the base level
        out = extrema_red(view(dts.input, i_start:i_stop))
    end
    return out
end

subselect_index(ndx, ndx_base) = ndx_base + ndx - 1
subselect_bin_bounds(bbounds::NTuple{2, <:Integer}, ndx_base::Integer) = subselect_index.(bbounds, ndx_base)

function downsamp_req(
    dts::CachingDynamicTs{S, A}, x_start, x_end, reqpts::Integer
) where {S, A}
    (in_range, i_begin, i_end) = downsamp_range_check(dts, x_start, x_end)
    if ! in_range
        return (Vector{Float64}(), Vector{NTuple{2,S}}(), false)
    end
    nbase = n_ndx(i_begin, i_end)
    binsize = max(fld(nbase, reqpts), 1)
    ncache = length(dts.cachearrays)
    decno = min(floor(Int, log10(binsize)), ncache)
    if  decno > 0
        was_downsampled = true
        nbin = cld(nbase, binsize)
        ys = Vector{NTuple{2, S}}(nbin)
        for binno in 1:nbin
            (bin_start, bin_stop) = subselect_bin_bounds(bin_bounds(binno, binsize, nbase), i_begin)
            ys[binno] = extrema_range(dts, bin_start, bin_stop)
        end
        bbounds = bin_bounds.(1:nbin, binsize, nbase)
    else
        was_downsampled = binsize > 1
        i_range = i_begin:i_end

        # Clip data
        subview = view(dts.input, i_range)
        mm = MaxMin(subview, binsize)
        ys = collect(mm)
        bbounds = bin_bounds(mm)
    end
    actual_start_x = ndx_to_t(i_begin, dts.fs, dts.offset)
    xs = ndx_to_t(bin_center(bbounds), dts.fs, actual_start_x)
    return (xs::Vector{Float64}, ys::Vector{NTuple{2,S}}, was_downsampled)
end

duration(dts::CachingDynamicTs) = duration(length(dts.input), dts.fs, dts.offset)

decade_ndx_conversion(i::Integer, dec::Integer) = fld(i - 1, 10 ^ dec) + 1
ndx_to_dec_ndx(args...) = decade_ndx_conversion(args...)
dec_ndx_to_ndx(i::Real, dec::Integer) = bin_center(i, 10 ^ dec)
dec_ndx_to_ndx(is::AbstractVector, dec::Integer) = bin_center.(is, 10 ^ dec)

struct MappedDynamicDownsampler{E, D<:DynamicDownsampler{E}} <: DynamicDownsampler{E}
    downsampler::D
    fmap::Function
end
function MappedDynamicDownsampler(d::D, fmap::Function) where {E, D<:DynamicDownsampler{E}}
    return MappedDynamicDownsampler{E,D}(d, fmap)
end

function downsamp_req(mdds::D, xb, xe, npts::Integer) where
    {D<:MappedDynamicDownsampler}
    (xs, ys, was_downsampled) = downsamp_req(mdds.downsampler, xb, xe, npts)
    mys = mdds.fmap(ys)
    return (xs, mys, was_downsampled)
end

duration(d::MappedDynamicDownsampler) = duration(d.downsampler)

make_shifter(shift) = (x) -> shift_extrema!(shift, x)

function shift_extrema(shift, ys::T) where {S, T<:NTuple{2,S}}
    (ys[1] + shift, ys[2] + shift)
end
shift_extrema(shift, y::Number) = y + shift

function shift_extrema!(shift, dest::B, ys::A) where
    {S, A<:AbstractVector{S}, B<:AbstractVector{S}}
    for (i, tup) in enumerate(ys)
        dest[i] = shift_extrema(shift, ys[i])
    end
end

function shift_extrema!(shift, ys)
    shift_extrema!(shift, ys, ys)
    return ys
end

function shift_extrema(shift, ys::AbstractVector)
    shifted = similar(ys)
    shift_extrema!(shift, shifted, ys)
    return shifted
end

