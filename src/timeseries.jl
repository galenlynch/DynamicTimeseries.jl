abstract type Downsampler{T, N} <: AbstractArray{T, N} end

immutable MaxMin{T<:Number, S, A<:AbstractArray} <: Downsampler{Tuple{T, T}, 1}
    input::A
    binsize::Int
    function MaxMin{T, S, A}(input::A, binsize::Integer) where {T<:Number, S<:Number, A<:AbstractVector{S}}
        return new(input, convert(Int, binsize))
    end
    function MaxMin{T, S, A}(input::A, binsize::Integer) where {T<:Number, S<:Number, A<:AbstractArray{S, 2}}
        @assert size(input, 1) == 2 "assumes max and min are on first dimension"
        return new(input, convert(Int, binsize))
    end
    function MaxMin{T, S, A}(input::A, binsize::Integer) where {T<:Number, S<:NTuple{2, T}, A<:AbstractVector{S}}
        return new(input, convert(Int, binsize))
    end
end
MaxMin(a::A, n::Integer) where {T<:Number, A<:AbstractVector{T}} = MaxMin{T, T, A}(a, n)
MaxMin(a::A, n::Integer) where {T<:Number, S<:NTuple{2, T}, A<:AbstractVector{S}} = MaxMin{T, S, A}(a, n)
MaxMin(a::A, n::Integer) where {T<:Number, A<:AbstractArray{T, 2}} = MaxMin{T, T, A}(a, n)

length(a::M) where {T, S, A<:AbstractVector, M<:MaxMin{T, S, A}} = cld(length(a.input), a.binsize)
length(a::M) where {T, S, A<:AbstractArray{T, 2}, M<:MaxMin{T, S, A}} = cld(size(a.input, 2), a.binsize)

size(a::MaxMin) = (length(a),)
binsize(a::MaxMin) = a.binsize

function getindex(a::M, i::Int) where {T<: Number, A<:AbstractVector, M<:MaxMin{T, T, A}}
    @boundscheck @assert checkbounds(Bool, a, i) "Index out of bounds"
    (idx_start, idx_stop) = bin_bounds(i, a.binsize, length(a.input))
    idx_range = idx_start:idx_stop
    return extrema(view(a.input, idx_range))
end
function getindex(a::M, i::Int) where {T<: Number, S<:NTuple{2, T}, A<:AbstractVector, M<:MaxMin{T, S, A}}
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
function getindex(a::M, i::Int) where {T<: Number, A<:AbstractArray{T, 2}, M<:MaxMin{T, T, A}}
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

setindex!(::MaxMin, ::Int) = throw(ReadOnlyMemoryError())

bin_bounds(a::MaxMin) = bin_bounds.(eachindex(a), a.binsize, length(a.input))
bin_bounds(i::Integer, a::MaxMin) = bin_bounds(i, a.binsize, length(a.input))

abstract type DynamicDownsampler end

immutable DynamicTs{T<:Number, A<:AbstractVector{T}} <: DynamicDownsampler
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
    offset::Real = 0
) where {T <: Number, A <: AbstractVector{T}}
    return DynamicTs{T, A}(input, convert(Float64, fs), convert(Float64, offset))
end

function downsamp_req(dts::DynamicTs, x_start::Real, x_end::Real, reqpoints::Integer)
    nin = length(dts.input)

    # Find bounding indices
    i_begin = clipind(x_to_ndx(x_start, dts.fs, dts.offset), nin)
    i_end = clipind(x_to_ndx(x_end, dts.fs, dts.offset), nin)

    # Calculate binsize
    npt = n_ndx(i_begin, i_end)
    binsize = convert(Int, ceil(npt / reqpoints))

    # Clip data
    subview = view(dts.input, i_begin:i_end)

    # Create downsampled data
    ys = MaxMin(subview, binsize)
    bin_start_x = ndx_to_x(i_begin, dts.fs, dts.offset)
    xs = ndx_to_x(bin_center(bin_bounds(ys)), dts.fs, bin_start_x)

    return (xs, ys)
end
function downsamp_req(dts::DynamicTs, x_start::Real, x_end::Real, reqpoints::AbstractFloat)
    return downsamp_req(dts, x_start, x_end, convert(Int, floor(reqpoints)))
end

immutable CachingDynamicTs{T<:Number, A<:AbstractVector{T}} <: DynamicDownsampler
    input::A
    fs::Float64
    offset::Float64
    cachefiles::Vector{Array{T, 2}}
    cachepaths::Vector{String}
    function CachingDynamicTs{T, A}(
        input::A,
        fs::Float64,
        offset::Float64,
        sizehint::Int,
        autoclean::Bool = true
    ) where {T<:Number, A<:AbstractVector{T}}
        nsamp = length(input)
        if sizehint < nsamp
            # Preallocate
            ndecade = convert(Int, ceil(log10(nsamp / sizehint)))
            cachefiles = Vector{Array{T, 2}}(ndecade)
            cachepaths = Vector{String}(ndecade)

            # Make first cache file
            (cachepaths[1], npair) = write_cache_file(input, autoclean)
            cachefiles[1] = open_cache_file(T, npair, cachepaths[1])
            # Make subsequent cache files
            for dno = 2:ndecade
                try
                    (cachepaths[dno], npair) = write_cache_file(cachefiles[dno - 1], autoclean)
                    cachefiles[dno] = open_cache_file(T, npair, cachepaths[dno])
                catch
                    for i = 1:(dno - 1)
                        rm(cachepaths[i])
                    end
                    rethrow()
                end
            end
        else
            cachefiles = Vector{Array{T, 2}}()
            cachepaths = Vector{String}()
        end
        return new(input, fs, offset, cachefiles, cachepaths)
    end
    function CachingDynamicTs{T, A}(
        input::A,
        fs::Float64,
        offset::Float64,
        cachepaths::Vector{String},
        cachelengths::Vector{Int},
        autoclean::Bool = true
    ) where {T<:Number, A<:AbstractVector{T}}
        # Preallocate
        ndecade = length(cachepaths)
        cachefiles = Vector{Array{T, 2}}(ndecade)
        try
            for dno in 1:ndecade
                cachefiles[dno] = open_cache_file(T, cachelengths[dno], cachepaths[dno])
                autoclean && atexit(() -> rm(cachepaths[dno]))
            end
        catch
            autoclean && foreach(rm, cachepaths)
            retrhow()
        end
        return new(input, fs, offset, cachefiles, cachepaths)
    end
end
function CachingDynamicTs(
    input::A,
    fs::Real,
    offset::Real = 0,
    sizehint::Integer = 300, # x dimension in pixels of a small window?
    autoclean::Bool = true
) where {T <: Number, A <: AbstractVector{T}}
    return CachingDynamicTs{T, A}(
        input,
        convert(Float64, fs),
        convert(Float64, offset),
        convert(Int, sizehint),
        autoclean
    )
end
function CachingDynamicTs(
    input::A,
    fs::Real,
    offset::Real,
    cachepaths::Vector{String},
    cachelengths::Vector{Int},
    autoclean::Bool = true
) where {T <: Number, A <: AbstractVector{T}}
    return CachingDynamicTs{T, A}(
        input,
        convert(Float64, fs),
        convert(Float64, offset),
        cachepaths,
        cachelengths,
        autoclean
    )
end

function write_cache_file(
    input::A,
    autoclean::Bool = true
) where {T<: Number, N, A<:AbstractArray{T, N}}
    mm = MaxMin(input, 10)
    (path, io) = mktemp()
    try
        for extremum in mm
            write(io, extremum...)
        end
    catch
        close(io)
        rm(path)
        rethrow()
    end
    close(io)
    autoclean && atexit(() -> rm(path))
    npair = length(mm)
    return (path, npair)
end

function open_cache_file(
    ::Type{T},
    npair::Integer,
    path::AbstractString
) where T
    cachearray = open(path, "r") do ior
        Mmap.mmap(ior, Array{T, 2}, (2, npair))
    end
    return cachearray
end

clean(a::CachingDynamicTs) = foreach(rm, a.cachepaths)

function scavenge_cache(a::CachingDynamicTs)
    ncache = length(a.cachefiles)
    cachelengths = Vector{Int}(ncache)
    for (i, cachefile) in enumerate(a.cachefiles)
        cachelengths[i] = size(cachefile, 2)
    end
    return (a.cachepaths, cachelengths)
end

function downsamp_req(
    dts::C,
    x_start::Real,
    x_end::Real,
    reqpts::Integer
) where {T<:Number, A<:AbstractVector, C<:CachingDynamicTs{T,A}}
    # Err on the side of too many
    nin = length(dts.input)
    i_begin = clipind(x_to_ndx(x_start, dts.fs, dts.offset), nin)
    i_end = clipind(x_to_ndx(x_end, dts.fs, dts.offset), nin)
    nbase = n_ndx(i_begin, i_end)
    ncache = length(dts.cachefiles)
    decno = min(convert(Int, floor(log10(nbase / reqpts))), ncache)
    if decno >= 1
        di_begin = decade_ndx_conversion(i_begin, decno)
        di_end = decade_ndx_conversion(i_end, decno)
        di_range = di_begin:di_end
        nsubsamp = n_ndx(di_begin, di_end)

        binsize = convert(Int, ceil(nsubsamp / reqpts))
        subview = view(dts.cachefiles[decno], di_range)
        ys = MaxMin(subview, binsize)

        dec_mult = 10 ^ decno
        x_decade_adj = ndx_to_x((1 + dec_mult) / 2, dts.fs)
        x_bin_adj = ndx_to_x(di_begin, dts.fs) * dec_mult
        bcs = bin_center(bin_bounds(ys)) * dec_mult
        xs = ndx_to_x(bcs, dts.fs, x_decade_adj + x_bin_adj + dts.offset)
    else
        binsize = convert(Int, ceil(nbase / reqpts))
        i_range = i_begin:i_end

        # Clip data
        subview = view(dts.input, i_range)
        ys = MaxMin(subview, binsize)

        bin_start_x = ndx_to_x(i_begin, dts.fs, dts.offset)
        xs = ndx_to_x(bin_center(bin_bounds(ys)), dts.fs, bin_start_x)
    end
    return (xs, ys)
end
function downsamp_req(dts::CachingDynamicTs, x_start::Real, x_end::Real, reqpoints::AbstractFloat)
    return downsamp_req(dts, x_start, x_end, convert(Int, floor(reqpoints)))
end

decade_ndx_conversion(i::Integer, dec::Integer) = fld(i - 1, 10 ^ dec) + 1
