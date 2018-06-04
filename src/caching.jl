const GL_CACHEPREFIX = "GLTSCACHE"

cache_reg(fid::Integer, ::Type{T}, n::Integer) where {T<:Number} = Regex(
    string(
        "^", GL_CACHEPREFIX, '_', fid, '_', T, '_', n, "_(\\d+)_(\\d+)_(\\d+)\$"
    )
)

function write_cache_file(
    ::Type{T},
    input::AbstractArray,
    autoclean::Bool = true,
    basename::AbstractString = tempname()
) where {T<:Downsampler}
    (cachelength, cachewidth, E, cachedata) = prepare_cachefile(T, input, basename)
    path = string(basename, '_', cachelength, '_', cachewidth)
    io = open(path, "w+")
    try
        write_cache_contents(io, cachedata)
    catch
        close(io)
        rm(path)
        rethrow()
    end
    close(io)
    autoclean && atexit(() -> rm(path))
    return (path, E, cachelength, cachewidth)
end

function write_cache_files(
    ::Type{D},
    input::A,
    sizehint::Integer,
    autoclean::Bool = true;
    fid::Integer = -1,
    cachedir::AbstractString = tempdir()
)where {D<:Downsampler, T<:Number, N, A<:AbstractArray{T, N}}
    nsamp = length(input)
    if fid > 0
        if ! isdir(cachedir)
            throw(ArgumentError(
                "cachedir ", cachedir, " is not a directory"
            ))
        end
        basestr = string(GL_CACHEPREFIX, '_', fid, '_', T, '_', 2) # 2 cache dims
        basename = joinpath(cachedir, basestr)
    else
        basename = string(tempname(), '_', T, '_', 2) # 2 cache dimensions
    end
    if sizehint < nsamp
        # Preallocate
        ndecade = convert(Int, ceil(log10(nsamp / sizehint)))
        cachepaths = Vector{String}(ndecade)
        lengths = Vector{Int}(ndecade)

        # Make first cache file
        dname = basename * "_1"
        (cachepaths[1], E, lengths[1], width) = write_cache_file(
            D, input, autoclean, dname
        )
        cachearr = open_cache_file(E, width, lengths[1], cachepaths[1])
        # Make subsequent cache files
        for dno = 2:ndecade
            dname = string(basename, '_', dno)
            try
                (cachepaths[dno], temp_el, lengths[dno], temp_width) = write_cache_file(
                    D, cachearr, autoclean, dname
                )
                temp_el != E && throw(Error("Element types are not all the same"))
                temp_width != width && throw(Error("Cache widths are not all the same"))
                cachearr = open_cache_file(E, width, lengths[dno], cachepaths[dno])
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
    return (cachepaths, E, lengths, width)
end

function open_cache_files(
    ::Type{E},
    width::Integer,
    lengths::AbstractArray{<:Integer},
    paths::Vector{String},
    autoclean::Bool = true
) where {E<:Number}
    ndecade = length(paths)
    cachearr = Vector{Array{E,2}}(ndecade)
    try
        for dno in 1:ndecade
            cachearr[dno] = open_cache_file(E, width, lengths[dno], paths[dno])
            autoclean && atexit(() -> rm(paths[dno]))
        end
    catch
        autoclean && foreach(rm, cachepaths)
        retrhow()
    end
    return cachearr
end

function open_cache_file(
    ::Type{E}, rows::Integer, cols::Integer, path::AbstractString
) where E<:Number
    open(path, "r") do ior
        Mmap.mmap(ior, Array{E, 2}, (rows, cols))
    end
end

function cacheinfo(a::CacheAccessor{<:Any,<:Any,<:Any,S}) where {S<:Number}
    cachelengths = map((x) -> size(x, 2), a.cachearrays)
    return (S, a.cachepaths, cachelengths)
end

function open_cache_files(
    ::Type{T}, cachedir::AbstractString, fid::Integer
) where {T<:Number}
    (fpaths, lengths, widths) = parse_cache_filenames(cachedir, fid, T, 2)
    if ! allsame(widths)
        throw(ErrorException("Widths must be the same"))
    end
    return open_cache_files(T, widths[1], lengths, fpaths, false)
end

function parse_cache_filenames(
    cachedir::AbstractString, fid::Integer, ::Type{T}, n::Integer
) where {T<:Number}
    matches = dir_match_files(cache_reg(fid, T, n), cachedir)
    nm = length(matches)

    lengths = Vector{Int}(nm)
    widths = Vector{Int}(nm)
    fpaths = Vector{String}(nm)
    if nm > 0
        filenos = Vector{Int}(nm)
        scratch = Vector{String}(nm)

        map!((x) -> x[1], scratch, matches)
        filenos .= parse.(Int, scratch)
        map!((x) -> x[2], scratch, matches)
        lengths .= parse.(Int, scratch)
        map!((x) -> x[3], scratch, matches)
        widths .= parse.(Int, scratch)
        map!((x) -> x.match, scratch, matches)
        fpaths .= joinpath.(cachedir, scratch)
        fpaths = fpaths[filenos]
        lengths = lengths[filenos]
    end
    return (fpaths, lengths, widths)
end

function validate_cache_arrays(
    cachepaths::AbstractVector{<:AbstractString},
    cachearrays::AbstractVector{<:AbstractArray},
    baselength::Integer,
    dim::Integer = 1,
    dec_factor::Integer = 10
)
    ncache = length(cachepaths)
    if ncache != length(cachearrays)
            throw(ArgumentError(
                "cachearrays and cachepaths must be the same length"
                ))
    end
    lens = Vector{Int}(ncache)
    lens .= size.(cachearrays, dim)
    last_len = baselength
    for l in lens
        if cld(last_len, l) != dec_factor
            throw(ArgumentError(
                "Cache arrays are not decreasing in size by a factor of $dec_factor"
            ))
        end
        last_len = l
    end
end

function sort_cache_files(
    cachepaths::AbstractArray{<:AbstractString},
    cachelengths::AbstractArray{<:Integer},
)
    p = sortperm(cachelengths; rev=true)
    sorted_cachelengths = cachelengths[p]
    sorted_cachepaths = cachepaths[p]
    return sorted_cachepaths, sorted_cachelengths
end
