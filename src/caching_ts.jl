struct CachingDynamicTs{
    W,E,T,C<:CacheAccessor{W,E,MaxMin,T}
} <: ExtremaDynDownsampler{T}
    accessor::C
end

function CachingDynamicTs(
    winput::W, cachearrays::Vector{Array{T, 2}}, cachepaths::Vector{String}
) where {T<:Number, W<:DynamicWindower{<:Any,T,1,<:Any}}
    println("Creating cache accessor")
    ca = CacheAccessor(winput, NTuple{2, T}, MaxMin, cachearrays, cachepaths)
    E = NTuple{2, T}
    C = CacheAccessor{W,E,MaxMin,T}
    CachingDynamicTs{W,E,T,C}(ca)
end

function CachingDynamicTs(
    winput::W,
    cachepaths::AbstractArray{<:AbstractString},
    cachelengths::AbstractArray{<:Integer},
    autoclean::Bool = false;
    checkfiles::Bool = true
) where {W<:DynamicWindower}
    if checkfiles
        println("Checking cachefiles")
        cachepaths, cachelengths = sort_cache_files(cachepaths, cachelengths)
    end
    println("Opening cache files")
    cachearrays = open_cache_files(
        cdts_dtype(W), cachepaths, cachelengths, autoclean
    )
    CachingDynamicTs(winput, cachearrays, cachepaths)
end

function CachingDynamicTs(
    winput::W, sizehint::Integer = 70, autoclean::Bool = true
) where {W<:DynamicWindower}
    println("Writing cache files")
    (cachepaths, cachelengths) = write_cache_files(
        cdts_dtype(W), basedata(winput), sizehint, autoclean
    )
    CachingDynamicTs(
        winput, cachepaths, cachelengths, false; checkfiles=false
    )
end

function CachingDynamicTs(
    input::A,
    fs::Real,
    offset::Real = 0,
    args...;
    kwargs...
) where {S<:Number, A<:AbstractVector{S}}
    CachingDynamicTs(DynamicWindower(input, fs, offset), args...; kwargs...)
end

function cdts_dtype(::Type{W}) where {T, W<:DynamicWindower{<:Any,T,1,<:Any}}
    E = NTuple{2, T}
    C = CacheAccessor{W,E,MaxMin,T}
    CachingDynamicTs{W, E, T, C}
end

cache_dims(::Type{D}) where D<:ExtremaDynDownsampler = 2

function open_cache_file(
    ::Type{D}, npair::Integer, path::AbstractString
) where {T<:Number, D<:ExtremaDynDownsampler{T}}
    cachearray = open(path, "r") do ior
        Mmap.mmap(ior, Array{T, 2}, (2, npair))
    end
    return cachearray
end

function new_cache_arrs(::Type{D}, n::Integer) where {T, D<:ExtremaDynDownsampler{T}}
    Vector{Array{T,2}}(n)
end

function prepare_cachefile(
    ::Type{T}, input::AbstractArray, basename::AbstractString
) where T<:ExtremaDynDownsampler
    mm = MaxMin(input, 10)
    npair = length(mm)
    return (npair, mm)
end

function write_cache_contents(::Type{T}, io::IO, mm::MaxMin) where T<:ExtremaDynDownsampler
    for extremum in mm
        write(io, extremum...)
    end
end

downsamp_req(dts::CachingDynamicTs, args...) = downsamp_req(dts.accessor, args...)

start_time(cdts::CachingDynamicTs) = start_time(cdts.accessor)
baselength(cdts::CachingDynamicTs) = baselength(cdts.accessor)
fs(cdts::CachingDynamicTs) = fs(cdts.accessor)
extrema(cdts::CachingDynamicTs) = downsamp_req(cdts, time_interval(cdts)..., 1)
