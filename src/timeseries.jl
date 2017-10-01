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

"Must implement downsamp_req and duration"
abstract type DynamicDownsampler end

immutable DynamicTs{A<:AbstractVector} <: DynamicDownsampler
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

function downsamp_req(dts::DynamicTs, x_start::Real, x_end::Real, reqpoints::Integer)
    nin = length(dts.input)

    # Find bounding indices
    i_begin = clipind(x_to_ndx(x_start, dts.fs, dts.offset), nin)
    i_end = clipind(x_to_ndx(x_end, dts.fs, dts.offset), nin)

    # Calculate binsize
    npt = n_ndx(i_begin, i_end)
    binsize = cld(npt, reqpoints)

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

duration(d::DynamicTs) = duration(length(d.input), d.fs, d.offset)

const GL_CACHEPREFIX = "GLTSCACHE"

immutable CachingDynamicTs{S<:Number, A<:AbstractVector{S}} <: DynamicDownsampler
    input::A
    fs::Float64
    offset::Float64
    cachefiles::Vector{Array{S, 2}}
    function CachingDynamicTs{S,A}(
        input::A,
        fs::Float64,
        offset::Float64,
        cachefiles::Vector{Array{S, 2}},
    ) where {S<:Number, A<:AbstractVector{S}}
        return new(input, fs, offset, cachefiles)
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

function write_cache_file(
    input::AbstractArray,
    autoclean::Bool = true,
    basename::AbstractString = tempname()
)
    mm = MaxMin(input, 10)
    npair = length(mm)
    path = basename * "_$npair"
    io = open(path, "w")
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
    return (path, npair)
end

function write_cache_files(
    input::A,
    sizehint::Integer,
    autoclean::Bool = true,
    basename::AbstractString = tempname()
)where {T<: Number, N, A<:AbstractArray{T, N}}
    nsamp = length(input)
    if sizehint < nsamp
        # Preallocate
        ndecade = convert(Int, ceil(log10(nsamp / sizehint)))
        cachepaths = Vector{String}(ndecade)
        lengths = Vector{Int}(ndecade)

        # Make first cache file
        dname = basename * "_1"
        (cachepaths[1], lengths[1]) = write_cache_file(input, autoclean, dname)
        cachearr = open_cache_file(T, lengths[1], cachepaths[1])
        # Make subsequent cache files
        for dno = 2:ndecade
            dname = basename * "_$dno"
            try
                (cachepaths[dno], lengths[dno]) = write_cache_file(cachearr, autoclean, dname)
                cachearr = open_cache_file(T, lengths[dno], cachepaths[dno])
            catch
                for i = 1:(dno - 1)
                    rm(cachepaths[i])
                end
                rethrow()
            end
        end
    else
        cachepaths = Vector{String}()
        lengths = Vector{Int}()
    end
    return (cachepaths, lengths)
end
function write_cache_files(cachedir::AbstractString, fid::Integer, autoclean::Bool = true, args...)
    basestr = string(GL_CACHEPREFIX, '_', fid)
    basename = joinpath(cachedir, basestr)
    write_cache_files(args..., autoclean, basename)
end

function open_cache_file(::Type{T}, npair::Integer, path::AbstractString) where T
    cachearray = open(path, "r") do ior
        Mmap.mmap(ior, Array{T, 2}, (2, npair))
    end
    return cachearray
end

function open_cache_files(
    ::Type{T},
    lengths::Vector{Int},
    paths::Vector{String},
    autoclean::Bool = true
) where T
    ndecade = length(paths)
    cachearr = Vector{Array{T, 2}}(ndecade)
    try
        for dno in 1:ndecade
            cachearr[dno] = open_cache_file(T, lengths[dno], paths[dno])
            autoclean && atexit(() -> rm(paths[dno]))
        end
    catch
        autoclean && foreach(rm, cachepaths)
        retrhow()
    end
    return cachearr
end

function open_cache_files(::Type{T}, cachedir::AbstractString, fid::Integer) where T
    fs = readdir()
    reg = Regex("^$(GL_CACHEPREFIX)_$(fid)_(\\d+)_(\\d+)\$")
    ms = match.(reg, fs)
    matches = filter((x) -> x != nothing, ms)
    nm = length(matches)

    filenos = Vector{Int}(nm)
    lengths = Vector{Int}(nm)
    scratch = Vector{String}(nm)
    fpaths = Vector{String}(nm)

    map!((x) -> x[1], scratch, matches)
    filenos .= parse.(Int, scratch)
    map!((x) -> x[2], scratch, matches)
    lengths .= parse.(Int, scratch)
    map!((x) -> x.match, scratch, matches)
    fpaths .= joinpath.(cachedir, scratch)
    fpaths = fpaths[filenos]
    lengths = lengths[filenos]
    return open_cache_files(T, lengths, fpaths, false)
end

function downsamp_req(dts::CachingDynamicTs, x_start::Real, x_end::Real, reqpts::Integer)
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

duration(dts::CachingDynamicTs) = duration(length(dts.input), dts.fs, dts.offset)

decade_ndx_conversion(i::Integer, dec::Integer) = fld(i - 1, 10 ^ dec) + 1

struct MappedDynamicDownsampler{D<:DynamicDownsampler} <: DynamicDownsampler
    downsampler::D
    fin::Function
    fx::Function
    fy::Function
end
function MappedDynamicDownsampler(d::DynamicDownsampler, fy::Function)
    MappedDynamicDownsampler(d, passthrough, passthrough, fy)
end

passthrough(args...) = args

function downsamp_req(
    mdds::D,
    xb::Real,
    xe::Real,
    npts::Integer
) where {D<:MappedDynamicDownsampler}
    trans_args = mdds.fin(xb, xe, npts)
    orig_res = downsamp_req(mdds.downsampler, trans_args...)
    println(orig_res)
    return (mdds.fx(orig_res...), mdds.fy(orig_res...))
end

duration(d::MappedDynamicDownsampler) = d.fx(duration(d.downsampler))
