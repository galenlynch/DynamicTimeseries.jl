const GL_CACHEPREFIX = "GLTSCACHE"

cache_reg(fid::Integer, ::Type{T}, n::Integer) where {T<:Number} = Regex(
    string("^", GL_CACHEPREFIX, '_', fid, '_', T, '_', n, "_(\\d+)_(\\d+)\$")
)

function write_cache_file(
    ::Type{T},
    input::AbstractArray,
    autoclean::Bool = true,
    basename::AbstractString = tempname()
) where {T<:ExtremaDynDownsampler}
    (cachelength, cachedata) = prepare_cachefile(T, input, basename)
    path = string(basename, '_', cachelength)
    io = open(path, "w+")
    try
        write_cache_contents(T, io, cachedata)
    catch
        close(io)
        rm(path)
        rethrow()
    end
    close(io)
    autoclean && atexit(() -> rm(path))
    return (path, cachelength)
end

function write_cache_files(
    ::Type{D},
    input::A,
    sizehint::Integer,
    autoclean::Bool = true;
    fid::Integer = -1,
    cachedir::AbstractString = tempdir()
)where {D<:DynamicDownsampler, T<:Number, N, A<:AbstractArray{T, N}}
    nsamp = length(input)
    cachedims = cache_dims(D)
    if fid > 0
        if ! isdir(cachedir)
            throw(ArgumentError(
                "cachedir ", cachedir, " is not a directory"
            ))
        end
        basestr = string(GL_CACHEPREFIX, '_', fid, '_', T, '_', cachedims)
        basename = joinpath(cachedir, basestr)
    else
        basename = string(tempname(), '_', T, '_', cachedims)
    end
    if sizehint < nsamp
        # Preallocate
        ndecade = convert(Int, ceil(log10(nsamp / sizehint)))
        cachepaths = Vector{String}(ndecade)
        lengths = Vector{Int}(ndecade)

        # Make first cache file
        dname = basename * "_1"
        (cachepaths[1], lengths[1]) = write_cache_file(D, input, autoclean, dname)
        cachearr = open_cache_file(D, lengths[1], cachepaths[1])
        # Make subsequent cache files
        for dno = 2:ndecade
            dname = string(basename, '_', dno)
            try
                (cachepaths[dno], lengths[dno]) = write_cache_file(
                    D, cachearr, autoclean, dname
                )
                cachearr = open_cache_file(D, lengths[dno], cachepaths[dno])
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

function open_cache_files(
    ::Type{D},
    paths::Vector{String},
    lengths::AbstractArray{<:Integer},
    autoclean::Bool = true
) where {D<:ExtremaDynDownsampler}
    ndecade = length(paths)
    cachearr = new_cache_arrs(D, ndecade)
    try
        for dno in 1:ndecade
            cachearr[dno] = open_cache_file(D, lengths[dno], paths[dno])
            autoclean && atexit(() -> rm(paths[dno]))
        end
    catch
        autoclean && foreach(rm, cachepaths)
        retrhow()
    end
    return cachearr
end

function cacheinfo(a::C) where {S<:Number, A, C<:CachingDynamicTs{S, A}}
    cachelengths = map((x) -> size(x, 2), a.cachearrays)
    return (S, a.cachepaths, cachelengths)
end

function open_cache_files(
    ::Type{D}, cachedir::AbstractString, fid::Integer
) where {T, D<:ExtremaDynDownsampler{T}}
    (fpaths, lengths) = parse_cache_filenames(cachedir, fid, T, 2)
    return open_cache_files(D, fpaths, lengths, false)
end

function parse_cache_filenames(
    cachedir::AbstractString, fid::Integer, ::Type{T}, n::Integer
) where {T<:Number}
    matches = dir_match_files(cache_reg(fid, T, n), cachedir)
    nm = length(matches)

    lengths = Vector{Int}(nm)
    fpaths = Vector{String}(nm)
    if nm > 0
        filenos = Vector{Int}(nm)
        scratch = Vector{String}(nm)

        map!((x) -> x[1], scratch, matches)
        filenos .= parse.(Int, scratch)
        map!((x) -> x[2], scratch, matches)
        lengths .= parse.(Int, scratch)
        map!((x) -> x.match, scratch, matches)
        fpaths .= joinpath.(cachedir, scratch)
        fpaths = fpaths[filenos]
        lengths = lengths[filenos]
    end
    return (fpaths, lengths)
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
