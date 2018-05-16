const GL_CACHEPREFIX = "GLTSCACHE"

cache_reg(fid::Integer) = Regex("^$(GL_CACHEPREFIX)_$(fid)_(\\d+)_(\\d+)\$")

function write_cache_file(
    input::AbstractArray,
    autoclean::Bool = true,
    basename::AbstractString = tempname()
)
    mm = MaxMin(input, 10)
    npair = length(mm)
    path = basename * "_$npair"
    io = open(path, "w+")
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
    return write_cache_files(args..., autoclean, basename)
end

function open_cache_file(::Type{T}, npair::Integer, path::AbstractString) where T
    cachearray = open(path, "r") do ior
        Mmap.mmap(ior, Array{T, 2}, (2, npair))
    end
    return cachearray
end

function open_cache_files(
    ::Type{T},
    paths::Vector{String},
    lengths::A,
    autoclean::Bool = true
) where {T, S<:Integer, A<:AbstractArray{S}}
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
    (fpaths, lengths) = parse_cache_filenames(cachedir, fid)
    return open_cache_files(T, fpaths, lengths, false)
end

function parse_cache_filenames(cachedir::AbstractString, fid::Integer)
    matches = dir_match_files(cache_reg(fid), cachedir)
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
