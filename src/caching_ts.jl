function prepare_cachefile(::Type{<:MaxMin}, input::AbstractArray, dim::Integer = 1)
    mm = MaxMin(input, 10)
    npair = length(mm)
    cachedims = (2, npair)
    return (cachedims::NTuple{2, Int}, mm)
end
