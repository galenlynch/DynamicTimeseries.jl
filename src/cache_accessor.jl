struct CacheAccessor{
    W<:DynamicWindower, E, D<:Downsampler, S<:Number, N
} <: AbstractDynamicDownsampler{E}
    winput::W
    cachearrays::Vector{Array{S,N}}
    cachepaths::Vector{String}

    function CacheAccessor{W,E,D,S,N}(
        winput::W,
        cachearrays::Vector{Array{S,N}},
        cachepaths::Vector{String}
    ) where {W<:DynamicWindower,E,D<:Downsampler,S<:Number,N}
        validate_cache_arrays(
            cachepaths,
            cachearrays,
            length(basedata(winput))
        )
        new(winput, cachearrays, cachepaths)
    end
end

function CacheAccessor(
    ::Type{D},
    winput::W,
    cachearrays::Vector{Array{S, N}},
    cachepaths::Vector{String},
    ::Type{E} = eltype_preview(D, A)
) where {
    D<:Downsampler,
    E,
    A,
    W<:DynamicWindower{<:Any, <:Any, <:Any, A},
    S<:Number,
    N
}
    CacheAccessor{W,E,D,S,N}(winput, cachearrays, cachepaths)
end

function CacheAccessor(
    ::Type{D},
    winput::DynamicWindower,
    ::Type{E},
    cachedims::AbstractVector{<:NTuple{N, <:Integer}},
    cachepaths::AbstractVector{<:AbstractString},
    autoclean::Bool = true,
    dim::Integer = N;
    checkfiles::Bool = true
) where {N, D<:Downsampler, E<:Number}
    if checkfiles
        cachepaths, cachedims = sort_cache_files(
            cachepaths, cachedims, dim
        )
    end
    cachearrays = open_cache_files(
        E, cachedims, cachepaths, autoclean
    )
    CacheAccessor(D, winput, cachearrays, cachepaths)
end

function CacheAccessor(
    ::Type{D},
    winput::DynamicWindower,
    sizehint::Integer = 70,
    autoclean::Bool = true
) where {D<:Downsampler}
    (cachepaths, E, cachedims) = write_cache_files(
        D, basedata(winput), sizehint, autoclean, winput.dim
    )
    CacheAccessor(
        D, winput, E, cachedims, cachepaths, false; checkfiles=false
    )
end

function CacheAccessor(
    ::Type{D},
    input::AbstractVector,
    fs::Real,
    args...;
    offset::Real = 0,
    dim::Integer = 1,
    f_overlap::Real = 0,
    wmin::Integer = 1,
    kwargs...
) where {D<:Downsampler}
    CacheAccessor(
        D,
        DynamicWindower(input, fs, offset, dim, f_overlap, wmin),
        args...;
        kwargs...
    )
end

function downsamp_req(
    dts::CacheAccessor{<:Any,E,D,<:Any},
    xb,
    xe,
    reqpts::Integer,
    exact::Bool = true,
    ::Type{T} = Int
) where {E, D<:Downsampler, T<:Integer}
    binsize, overlap, i_begin, i_end = downsamp_req_window_info(
        dts.winput, xb, xe, T(reqpts)
    )
    nbase = n_ndx(i_begin, i_end)

    if nbase <= 0 || reqpts == 0
        return (Vector{Float64}(), Vector{E}(), false)
    end

    ncache = T(length(dts.cachearrays))
    decno = min(floor(T, log10(binsize)), ncache)
    x_start = ndx_to_t(i_begin, fs(dts.winput), start_time(dts.winput))
    if  decno > 0
        was_downsampled = true
        if exact
            times, ys = exact_downsample(
                dts, nbase, binsize, i_begin, i_end, x_start
            )
        else
            times, ys = approximate_downsample(dts, binsize, i_begin, i_end, decno)
        end
    else
        times, wa, was_downsampled = make_windowedarray(
            dts.winput, binsize, overlap, i_begin, i_end
        )
        downsampler = D(wa)
        ys = collect(downsampler)
    end
    return (times, ys, was_downsampled)
end

start_time(ca::CacheAccessor) = start_time(ca.winput)
fs(ca::CacheAccessor) = fs(ca.winput)
baselength(ca::CacheAccessor) = baselength(ca.winput)

function approximate_downsample(
    dts::CacheAccessor{<:Any, <:Any, D, <:Any, N},
    binsize::T,
    i_begin::T,
    i_end::T,
    decno::T
) where {T<:Integer, D, N}
    di_begin = ndx_to_dec_ndx(i_begin, decno)
    di_end = ndx_to_dec_ndx(i_end, decno)
    nsubsamp = n_ndx(di_begin, di_end)

    slice_idx = make_slice_idx(N, N, di_begin:di_end)
    cache_slice = view(dts.cachearrays[decno], slice_idx...)

    dec_binsize = fld(binsize, T(10) ^ decno)
    wa = WindowedArray(cache_slice, dec_binsize, N)

    bbounds = bin_bounds(wa, T)
    for (i, bpair) in enumerate(bbounds)
        bbounds[i] = dec_ndx_to_ndx(bpair .+ di_begin .- one(T), decno)
    end
    bcenters = bin_center(bbounds)
    xs = ndx_to_t(bcenters, fs(dts.winput), start_time(dts.winput))

    ys = collect(D(wa))
    return (xs, ys)
end

function exact_downsample(
    dts::CacheAccessor{<:Any,E,<:Any,<:Any,<:Any},
    nbase::Integer,
    binsize::Integer,
    i_begin::Integer,
    i_end::Integer,
    x_start::Real
) where {E}
    nbin = cld(nbase, binsize)
    ys = Vector{E}(nbin)
    for binno in 1:nbin
        (bin_start, bin_stop) = subselect_bin_bounds(
            bin_bounds(binno, binsize, nbase), i_begin
        )
        ys[binno], _ = reduce_downsample_caches(dts, bin_start, bin_stop)
    end
    bbounds = Vector{NTuple{2, Int}}(nbin)
    bbounds .= bin_bounds.(1:nbin, binsize, nbase)
    xs = ndx_to_t(bin_center(bbounds), fs(dts.winput), x_start)
    return (xs, ys)
end

"Recursively find downsampled values using cached arrays"
function reduce_downsample_caches(
    dts::CacheAccessor{W,E,D,<:Any,N},
    i_start::Integer,
    i_stop::Integer
) where {M, W<:DynamicWindower{<:Any,<:Any,M,<:Any}, E, D<:Downsampler, N}
    # In order to find the downsampled value that you would get without taking
    # advantage of the cached downsampled values, we can combine cached values
    # together and traverse the caches (as well as the input data) to get the
    # full answer.

    nbase = n_ndx(i_start, i_stop)
    decno = min(floor(Int, log10(nbase)), length(dts.cachearrays))
    D_concrete = concrete_type(D, windowtype(W))
    if decno > 0
        cached_ds = Vector{E}(3) # storage for left, right, and center
        weights = Vector{Int}(3)
        n_input = length(basedata(dts.winput)) # number of original samples


        # Find the indices of cached values that are completely contained
        # in this slice
        contained_first = dec_ndx_greater(i_start, decno)
        if n_input == i_stop # last index
            contained_last = size(dts.cachearrays[decno], N)
        else
            contained_last = dec_ndx_lesser(i_stop, decno)
        end

        slice_idx = make_slice_idx(N, N, contained_first:contained_last)
        contained_view = view(dts.cachearrays[decno], slice_idx...)

        # Reduce contained cache values
        if size(contained_view, N) > 0 # if the view is not empty
            # Reduce downsampled values from this cache level
            reduce_weights = fill(10 ^ decno, size(contained_view, N))
            cached_ds[1], weights[1] = downsamp_reduce_cache(
                D_concrete, contained_view, reduce_weights
            )
            extrema_i = 2 # advance cached_ds index
            remainder_left_stop = dec_ndx_2_base_start(contained_first, decno) - 1
            remainder_right_start = dec_ndx_2_base_end(contained_last, decno) + 1
        else
            # This cache level will not be used, skip it
            extrema_i = 1 # Do not advance the index
            deleteat!(cached_ds, 1) # Shrink cahced_ds by one
            deleteat!(weights, 1)
            # Split indicies in a cache-friendly manner
            next_dec = decno - 1
            remainder_left_stop = (10 ^ next_dec) * ndx_to_dec_ndx(
                i_start + fld(nbase, 2),
                next_dec
            ) # Should be at the end of a cache
            remainder_right_start = remainder_left_stop + 1
        end

        # Recurse on 'left' remainder
        if i_start <= remainder_left_stop
            cached_ds[extrema_i], weights[extrema_i] = reduce_downsample_caches(
                dts, i_start, remainder_left_stop
            )
            extrema_i += 1
        else
            deleteat!(cached_ds, extrema_i)
            deleteat!(weights, extrema_i)
        end

        # Recurse on 'right' remainder
        if i_stop >= remainder_right_start
            cached_ds[extrema_i], weights[extrema_i] = reduce_downsample_caches(
                dts, remainder_right_start, i_stop
            )
            extrema_i += 1
        else
            deleteat!(cached_ds, extrema_i)
            deleteat!(weights, extrema_i)
        end

        out, weight = downsamp_reduce_cache(D_concrete, cached_ds, weights)
    else # We're at the base level
        slice_idx = make_slice_idx(M, dts.winput.dim, i_start:i_stop)
        baseview = view(basedata(dts.winput), slice_idx...)
        reduce_weights = fill(1, size(baseview, dts.winput.dim))
        out, weight = downsamp_reduce(
            D_concrete, baseview, reduce_weights, dts.winput.dim
        )
    end
    return out, weight
end

function write_cache_contents(io::IO, ds::Downsampler)
    for downsampled in ds
        write(io, downsampled...)
    end
end

subselect_index(ndx, ndx_base) = ndx_base + ndx - 1
function subselect_bin_bounds(bbounds::NTuple{2, <:Integer}, ndx_base::Integer)
    subselect_index.(bbounds, ndx_base)
end

dec_ndx_greater(i, dec) = cld(i - 1, 10 ^ dec) + 1
dec_ndx_lesser(i, dec) = fld(i, 10 ^ dec)
dec_ndx_2_base_start(i, dec) = (i - 1) * 10 ^ dec + 1
dec_ndx_2_base_end(i, dec) = i * 10 ^ dec

function ndx_to_dec_ndx(i::T, dec::T) where {T<:Integer}
    fld(i - one(T), T(10) ^ dec) + one(T)
end

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
