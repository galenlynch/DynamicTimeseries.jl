function prepare_cachefile(::Type{<:MaxMin}, input::AbstractArray)
    mm = MaxMin(input, 10)
    npair = length(mm)
    cachedims = (2, npair)
    return (cachedims::NTuple{2,Int}, mm)
end

function prepare_next_cachefile(::Type{M}, args...) where {M<:MaxMin}
    prepare_cachefile(M, args...)
end

function extrema(ca::CacheAccessor{<:Any,<:Any,<:MaxMin,<:Any,<:Any})
    _, ys, _ = downsamp_req(ca, time_interval(ca)..., 1, false)
    extrema_red(ys)
end
