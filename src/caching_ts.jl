function prepare_cachefile(
    ::Type{D}, input::A, dim::Integer = 1
) where {D<:MaxMin, E<:Number, A<:AbstractArray{E, <:Any}}
    println("Size of input is ", size(input), " and dim is ", dim)
    mm = MaxMin(input, 10)
    npair = length(mm)
    cachedims = (2, npair)
    println("returning cachedims ", cachedims)
    return (cachedims::NTuple{2, Int}, E, mm)
end
