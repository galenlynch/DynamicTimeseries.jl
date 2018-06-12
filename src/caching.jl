const GL_CACHEPREFIX = "GLTSCACHE"

cache_reg(fid::Integer, ::Type{T}, n::Integer) where {T<:Number} = Regex(
    string(
        "^", GL_CACHEPREFIX, '_', fid, '_', T, '_', "_(\\d+)_([\\d_]+)\$"
    )
)

function write_cache_file(
    ::Type{D},
    input::AbstractArray{<:Any, N},
    autoclean::Bool = true,
    basename::AbstractString = tempname(),
    dim::Integer = N
) where {N, D<:Downsampler}
    (cachedims, cachedata) = prepare_cachefile(D, input, dim)
    ndim = length(cachedims)
    dim_els = Vector{String}(ndim + 1)
    dim_els[1] = string(ndim)
    for i in 1:ndim
        dim_els[i + 1] = string(cachedims[i])
    end
    dim_str = join(dim_els, '_')
    path = string(basename, '_', dim_str)
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
    return (path, cachedims)
end

function write_cache_files(
    ::Type{D},
    input::A,
    sizehint::Integer,
    autoclean::Bool = true,
    dim::Integer = 1;
    fid::Integer = -1,
    cachedir::AbstractString = tempdir()
)where {D<:Downsampler, A<:AbstractArray}
    input_eltype = arr_eltype_preview(D, A)
    nsamp = length(input)
    if fid > 0
        if ! isdir(cachedir)
            throw(ArgumentError(
                "cachedir ", cachedir, " is not a directory"
            ))
        end
        basestr = string(GL_CACHEPREFIX, '_', fid, '_', input_eltype)
        basename = joinpath(cachedir, basestr)
    else
        basename = string(tempname(), '_', input_eltype)
    end
    E = arr_eltype_preview(D, A)
    if sizehint < nsamp
        # Preallocate
        ndecade = convert(Int, ceil(log10(nsamp / sizehint)))
        cachepaths = Vector{String}(ndecade)

        # Make first cache file
        dname = basename * "_1"
        (cachepaths[1], cachedim) = write_cache_file(
            D, input, autoclean, dname, dim
        )
        cachedims = Vector{typeof(cachedim)}(ndecade)
        cachedims[1] = cachedim
        cachearr = open_cache_file(E, cachedims[1], cachepaths[1])
        # Make subsequent cache files
        for dno = 2:ndecade
            dname = string(basename, '_', dno)
            try
                (cachepaths[dno], cachedims[dno]) = write_cache_file(
                    D, cachearr, autoclean, dname
                )
                cachearr = open_cache_file(E, cachedims[dno], cachepaths[dno])
            catch
                for i = 1:(dno - 1)
                    rm(cachepaths[i])
                end
                rethrow()
            end
        end
    else
        cachepaths = Vector{String}()
        cachedims = Vector{NTuple{2, Int}}()
    end
    return (cachepaths, E, cachedims)
end

arr_bit_eltype(::Type{A}) where {T<:Number, A<:AbstractArray{T, <:Any}} = T

function arr_bit_eltype(
    ::Type{B}
) where {T<:Number, A<:AbstractArray{T, <:Any}, B<:AbstractVector{A}}
    T
end

function write_cache_files(
    ::Type{D},
    winput::DynamicWindower,
    sizehint::Integer,
    autoclean::Bool = true;
    fid::Integer = -1,
    cachedir::AbstractString = tempdir()
) where D <: Downsampler
    write_cache_files(
        D, basedata(winput), sizehint, autoclean, winput.dim;
        fid = fid, cachedir = cachedir
    )
end

function open_cache_files(
    ::Type{E},
    dims::AbstractVector{<:NTuple{N, <:Integer}},
    paths::AbstractVector{<:AbstractString},
    autoclean::Bool = true
) where {N, E<:Number}
    ndecade = length(paths)
    length(dims) == ndecade || error("paths and dims must be same length")
    cachearr = Vector{Array{E,N}}(ndecade)
    try
        for dno in 1:ndecade
            cachearr[dno] = open_cache_file(E, dims[dno], paths[dno])
            autoclean && atexit(() -> rm(paths[dno]))
        end
    catch
        autoclean && foreach(rm, cachepaths)
        retrhow()
    end
    return cachearr
end

function open_cache_files(
    ::Type{T}, cachedir::AbstractString, fid::Integer
) where {T<:Number}
    (fpaths, dims) = parse_cache_filenames(cachedir, fid, T, 2)
    return open_cache_files(T, dims, fpaths, false)
end

function open_cache_file(
    ::Type{E}, dims::NTuple{N, <:Integer}, path::AbstractString
) where {N, E<:Number}
    open(path, "r") do ior
        Mmap.mmap(ior, Array{E, N}, dims)
    end
end

function cacheinfo(a::CacheAccessor{<:Any,<:Any,<:Any,S}) where {S<:Number}
    cachedims = map(size, a.cachearrays)
    return (S, a.cachepaths, cachedims)
end

function parse_cache_filenames(
    cachedir::AbstractString, fid::Integer, ::Type{T}, n::Integer
) where {T<:Number}
    matches = dir_match_files(cache_reg(fid, T, n), cachedir)
    nm = length(matches)

    dims = Vector{NTuple{n, Int}}(nm)
    dim_scratch = Vector{Int}(n)
    fpaths = Vector{String}(nm)
    if nm > 0
        filenos = Vector{Int}(nm)
        for (i, m) in enumerate(matches)
            filenos[i] = parse(Int, m[1])
            dim_strs = split(m[2], '_')
            length(dim_strs) == n || error("Did not find the correct number of dims")
            for j in 1:n
                dim_scratch[j] = parse(Int, dim_strs[j])
            end
            dims[i] = (dim_scratch...)
            fpaths[i] = joinpath(cachedir, m.match)
        end
        p = sortperm(filenos)
        fpaths = fpaths[p]
        dims = dims[p]
    end
    return (fpaths, dims)
end

function validate_cache_arrays(
    cachepaths::AbstractVector{<:AbstractString},
    cachearrays::AbstractVector{<:AbstractArray{<:Any, N}},
    baselength::Integer,
    dim::Integer = N,
    dec_factor::Integer = 10
) where N
    ncache = length(cachepaths)
    if ncache != length(cachearrays)
            throw(ArgumentError(
                "cachearrays and cachepaths must be the same length"
                ))
    end
    last_len = baselength
    for carr in cachearrays
        l = size(carr, dim)
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
    cachedims::AbstractArray{N, <:Integer},
    dim::Integer
) where {N}
    cachlengths = size.(cachedims, dim)
    p = sortperm(cachelengths; rev=true)
    sorted_cachedims = cachedims[p]
    sorted_cachepaths = cachepaths[p]
    return sorted_cachepaths, sorted_cachedims
end
