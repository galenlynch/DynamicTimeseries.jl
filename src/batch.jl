function downsamp_batch_mmap(
    ::Type{DownsamplerType}, # downsampler type
    arrs::AbstractVector{<:AbstractArray},
    paths::AbstractVector{<:AbstractString},
    fss::Union{<:Real,AbstractVector{<:Real}},
    offbytes::Union{<:Integer,AbstractVector{<:Integer}} = 0,
    t_offsets::AbstractVector{<:Real} = Int[],
    sizehint::Integer = 70;
    cachedir::AbstractString = tempdir(),
    fids::AbstractVector{<:Integer} = Int[],
    autoclean::Bool = true,
) where {DownsamplerType<:Downsampler}
    t_offs = isempty(t_offsets) ? zeros(Int, length(arrs)) : t_offsets
    samplerates = isa(fss, Real) ? fill(fss, length(arrs)) : fss
    mm_offs = isa(offbytes, Integer) ? fill(offbytes, length(arrs)) : offbytes
    mmap_types = typeof.(arrs)
    mmap_sizes = size.(arrs)
    cpaths, ctypes, cdims = cachefiles_batch_mmap(
        DownsamplerType,
        paths,
        mmap_types,
        mmap_sizes,
        mm_offs,
        sizehint,
        cachedir,
        fids,
        false,
    )
    CacheAccessor.(
        DownsamplerType,
        arrs,
        samplerates,
        t_offs,
        ctypes,
        cdims,
        cpaths,
        autoclean,
    )
end

function downsamp_batch_mmap(
    ::Type{DownsamplerType}, # downsampler type
    arrs::AbstractVector{<:AbstractArray},
    fss::Union{<:Real,AbstractVector{<:Real}},
    args...,
) where {DownsamplerType<:Downsampler}
    na = length(arrs)
    na > 0 || throw(ArgumentError("arrs cannot be empty"))

    outs = to_mmap.(arrs)

    mm_arrs = Vector{typeof(outs[1][1])}(undef, na)
    offbytes = zeros(Int, na)
    paths = Vector{String}(undef, na)

    for (i, (a, p)) in enumerate(outs)
        mm_arrs[i] = a
        paths[i] = p
    end

    downsamp_batch_mmap(DownsamplerType, mm_arrs, paths, fss, offbytes, args...)
end

function cachefiles_batch_mmap(
    ::Type{DownsamplerType},
    paths::AbstractVector{<:AbstractString},
    mmap_types::AbstractVector{DataType},
    mmap_sizes::AbstractVector{<:NTuple{<:Any,<:Integer}},
    offbs::Union{Integer,AbstractVector{<:Integer}},
    sizehint::Integer,
    cachedir::AbstractString = tempdir(),
    fids::AbstractArray{<:Integer} = Int[],
    autoclean::Bool = true,
) where {DownsamplerType<:Downsampler}
    np = length(paths)
    np > 0 || throw(ArgumentError("Paths cannot be empty"))
    if isempty(fids)
        fids = fill(-1, np)
    elseif length(fids) != np
        throw(ArgumentError("fids must be empty or the same size as paths"))
    end
    mm_offs = isa(offbs, Integer) ? fill(offbs, np) : offbs

    outs = pmap(
        (p, mt, ms, ob, fid) -> cachefiles_mmap_remote(
            DownsamplerType,
            p,
            mt,
            ms,
            ob,
            sizehint,
            cachedir,
            fid,
        ),
        paths,
        mmap_types,
        mmap_sizes,
        mm_offs,
        fids,
    )

    # Destructure pmap output
    cpaths = Vector{Vector{String}}(undef, np)
    ctypes = Vector{DataType}(undef, np)
    csizes = Vector{typeof(outs[1][3])}(undef, np)
    for (i, (cp, ct, cs)) in enumerate(outs)
        cpaths[i] = cp
        autoclean && atexit(() -> rm.(cp))
        ctypes[i] = ct
        csizes[i] = cs
    end

    cpaths, ctypes, csizes
end

function cachefiles_mmap_remote(
    ::Type{DownsamplerType},
    path::AbstractString,
    ::Type{MmapType},
    mmap_size::NTuple{<:Any,<:Integer},
    offbytes::Integer,
    sizehint::Integer,
    cachedir::AbstractString = tempdir(),
    fid::Integer = -1,
) where {DownsamplerType<:Downsampler,MmapType<:AbstractArray}
    arr = open(path, "r") do io
        Mmap.mmap(io, MmapType, mmap_size, offbytes)
    end
    write_cache_files(DownsamplerType, arr, sizehint, false; cachedir = cachedir, fid = fid)
end
