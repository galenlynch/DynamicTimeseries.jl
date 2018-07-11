function downsamp_batch_mmap(
    ::Type{DownsamplerType}, # downsampler type
    arrs::AbstractVector{<:AbstractArray},
    paths::AbstractVector{<:AbstractString},
    fss::Union{<:Real, AbstractVector{<:Real}},
    offbytes::Union{<:Integer, AbstractVector{<:Integer}} = 0,
    t_offsets::AbstractVector{<:Real} = Int[],
    sizehint::Integer = 70
) where {DownsamplerType<:Downsampler}
    t_offs = isempty(t_offsets) ? zeros(Int, length(arrs)) : t_offsets
    samplerates = isa(fss, Real) ? fill(fss, length(arrs)) : fss
    mm_offs = isa(offbytes, Integer) ? fill(offbytes, length(arrs)) : offbytes
    mmap_types = typeof.(arrs)
    mmap_sizes = size.(arrs)
    cpaths, ctypes, cdims = cachefiles_batch_mmap(
        DownsamplerType, paths, mmap_types, mmap_sizes, mm_offs, sizehint
    )
    CacheAccessor.(
        DownsamplerType, arrs, samplerates, t_offs, ctypes, cdims, cpaths, true
    )
end
function downsamp_batch_mmap(
    ::Type{DownsamplerType}, # downsampler type
    arrs::AbstractVector{<:AbstractArray},
    fss::Union{<:Real, AbstractVector{<:Real}},
    args...
) where {DownsamplerType <: Downsampler}
    na = length(arrs)
    na > 0 || throw(ArgumentError("arrs cannot be empty"))

    outs = to_mmap.(arrs)

    mm_arrs = Vector{typeof(outs[1][1])}(na)
    offbytes = zeros(Int, na)
    paths = Vector{String}(na)

    for (i, (a, p)) in enumerate(outs)
        mm_arrs[i] = a
        paths[i] = p
    end

    downsamp_batch_mmap(
        DownsamplerType,
        mm_arrs,
        paths,
        fss,
        offbytes,
        args...
    )
end

function cachefiles_batch_mmap(
    ::Type{DownsamplerType},
    paths::AbstractVector{<:AbstractString},
    mmap_types::AbstractVector{DataType},
    mmap_sizes::AbstractVector{<:NTuple{<:Any, <:Integer}},
    offbs::AbstractVector{<:Integer},
    sizehint::Integer
) where {DownsamplerType<:Downsampler}
    np = length(paths)
    np > 0 || throw(ArgumentError("Paths cannot be empty"))

    outs = pmap(
        (p, mt, ms, ob) -> cachefiles_mmap_remote(
            DownsamplerType, p, mt, ms, ob, sizehint
        ),
        paths, mmap_types, mmap_sizes, offbs
    )

    # Destructure pmap output
    cpaths = Vector{Vector{String}}(np)
    ctypes = Vector{DataType}(np)
    csizes = Vector{typeof(outs[1][3])}(np)
    for (i, (cp, ct, cs)) in enumerate(outs)
        cpaths[i] = cp
        ctypes[i] = ct
        csizes[i] = cs
    end

    cpaths, ctypes, csizes
end

function cachefiles_mmap_remote(
    ::Type{DownsamplerType},
    path::AbstractString,
    ::Type{MmapType},
    mmap_size::NTuple{<:Any, <:Integer},
    offbytes::Integer,
    sizehint::Integer
) where {DownsamplerType <: Downsampler, MmapType <: AbstractArray}
    arr = open(path, "r") do io
        Mmap.mmap(io, MmapType, mmap_size, offbytes)
    end
    write_cache_files(DownsamplerType, arr, sizehint)
end
