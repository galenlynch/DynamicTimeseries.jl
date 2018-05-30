struct CachingDynamicTs{S<:Number, A<:AbstractVector{S}} <: ExtremaDownsampler{S}
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
        validate_cache_arrays(cachepaths, cachearrays, length(input), 2)
        new(input, fs, offset, cachearrays, cachepaths)
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
        cachepaths, cachelengths = sort_cache_files(cachepaths, cachelengths)
    end
    cachearrs = open_cache_files(
        CachingDynamicTs{S,A}, cachepaths, cachelengths, autoclean
    )
    CachingDynamicTs(input, fs, offset, cachearrs, cachepaths)
end

function CachingDynamicTs(
    input::A,
    fs::Real,
    offset::Real = 0,
    sizehint::Integer = 70, # x dimension in pixels of a small window?
    autoclean::Bool = true
) where {S<:Number, A<:AbstractArray{S}}
    (cachepaths, cachelengths) = write_cache_files(
        CachingDynamicTs{S,A}, input, sizehint, autoclean
    )
    CachingDynamicTs(
        input, fs, offset, cachepaths, cachelengths, false; checkfiles=false
    )
end

cache_dims(::Type{D}) where D<:ExtremaDownsampler = 2

function open_cache_file(
    ::Type{D}, npair::Integer, path::AbstractString
) where {T<:Number, D<:ExtremaDownsampler{T}}
    cachearray = open(path, "r") do ior
        Mmap.mmap(ior, Array{T, 2}, (2, npair))
    end
    return cachearray
end

function new_cache_arrs(::Type{D}, n::Integer) where {T, D<:ExtremaDownsampler{T}}
    Vector{Array{T,2}}(n)
end

function prepare_cachefile(
    ::Type{T}, input::AbstractArray, basename::AbstractString
) where T<:ExtremaDownsampler
    mm = MaxMin(input, 10)
    npair = length(mm)
    return (npair, mm)
end

function write_cache_contents(::Type{T}, io::IO, mm::MaxMin) where T<:ExtremaDownsampler
    for extremum in mm
        write(io, extremum...)
    end
end

dec_ndx_greater(i, dec) = cld(i - 1, 10 ^ dec) + 1
dec_ndx_lesser(i, dec) = fld(i, 10 ^ dec)
dec_ndx_2_base_start(i, dec) = (i - 1) * 10 ^ dec + 1
dec_ndx_2_base_end(i, dec) = i * 10 ^ dec

"Recursively find extrema using cached arrays"
function extrema_range(
    dts::CachingDynamicTs{S, A}, i_start::Integer, i_stop::Integer
) where {S, A}
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
            remainder_left_stop = (10 ^ next_dec) * ndx_to_dec_ndx(
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
function subselect_bin_bounds(bbounds::NTuple{2, <:Integer}, ndx_base::Integer)
    subselect_index.(bbounds, ndx_base)
end

function downsamp_req(
    dts::CachingDynamicTs{S, A},
    x_start_req,
    x_end_req,
    reqpts::Integer,
    exact::Bool=true,
    ::Type{T} = Int32
) where {S, A, T<:Integer}
    (in_range, i_begin, i_end) = downsamp_range_check(
        dts, x_start_req, x_end_req, T
    )
    nbase = n_ndx(i_begin, i_end)

    if ! in_range || nbase == 0 || reqpts == 0
        return (Vector{Float64}(), Vector{NTuple{2,S}}(), false)
    end

    binsize = max(fld(nbase, convert(T, reqpts)), one(T))
    ncache = T(length(dts.cachearrays))
    decno = min(floor(T, log10(binsize)), ncache)
    x_start = ndx_to_t(i_begin, dts.fs, dts.offset)
    if  decno > 0
        was_downsampled = true
        if exact
            xs, ys = exact_downsample(dts, nbase, binsize, i_begin, i_end, x_start)
        else
            xs, ys = approximate_downsample(dts, binsize, i_begin, i_end, decno)
        end
    else
        was_downsampled = binsize > 1
        i_range = i_begin:i_end

        # Clip data
        subview = view(dts.input, i_range)
        mm = MaxMin(subview, binsize)
        ys = collect(mm)
        bbounds = bin_bounds(mm)
        xs = ndx_to_t(bin_center(bbounds), dts.fs, x_start)
    end
    return (xs::Vector{Float64}, ys::Vector{NTuple{2,S}}, was_downsampled::Bool)
end

function exact_downsample(
    dts::CachingDynamicTs{S, <:Any},
     nbase::Integer,
     binsize::Integer,
     i_begin::Integer,
     i_end::Integer,
     x_start::Real
) where {S}
    nbin = cld(nbase, binsize)
    ys = Vector{NTuple{2, S}}(nbin)
    for binno in 1:nbin
        (bin_start, bin_stop) = subselect_bin_bounds(
            bin_bounds(binno, binsize, nbase), i_begin
        )
        ys[binno] = extrema_range(dts, bin_start, bin_stop)
    end
    bbounds = Vector{NTuple{2, Int}}(nbin)
    bbounds .= bin_bounds.(1:nbin, binsize, nbase)
    xs = ndx_to_t(bin_center(bbounds), dts.fs, x_start)
    return (xs, ys)
end

function approximate_downsample(
    dts::CachingDynamicTs,
    binsize::T,
    i_begin::T,
    i_end::T,
    decno::T
) where {T<:Integer}
    di_begin = ndx_to_dec_ndx(i_begin, decno)
    di_end = ndx_to_dec_ndx(i_end, decno)
    nsubsamp = n_ndx(di_begin, di_end)

    subview = view(dts.cachearrays[decno], 1:2, di_begin:di_end)

    dec_binsize = fld(binsize, T(10) ^ decno)
    mm = MaxMin(subview, dec_binsize)

    bbounds = bin_bounds(mm, T)
    for (i, bpair) in enumerate(bbounds)
        bbounds[i] = dec_ndx_to_ndx(bpair .+ di_begin .- one(T), decno)
    end
    bcenters = bin_center(bbounds)
    xs = ndx_to_t(bcenters, dts.fs, dts.offset)
    ys = collect(mm)
    return (xs, ys)
end

fs(dts::CachingDynamicTs) = dts.fs
start_time(dts::CachingDynamicTs) = dts.offset
baselength(dts::CachingDynamicTs) = length(dts.input)

ndx_to_dec_ndx(i::T, dec::T) where {T<:Integer} = fld(i - one(T), T(10) ^ dec) + one(T)

dec_ndx_to_ndx(i::T, dec::T) where {T<:Integer} = bin_bounds(i, T(10) ^ dec)
function dec_ndx_to_ndx(is::NTuple{2, <:Integer}, dec::Integer)
    (start_start, start_stop) = dec_ndx_to_ndx(is[1], dec)
    (stop_start, stop_stop) = dec_ndx_to_ndx(is[2], dec)
    return (start_start, stop_stop)
end
function dec_ndx_to_ndx(is::AbstractVector{T}, dec::T) where T<:Integer
    rs = Vector{NTuple{2, T}}(length(is))
    @. rs = dec_ndx_to_ndx(is, T(10) ^ dec)
end
