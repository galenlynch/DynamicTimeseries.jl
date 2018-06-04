function prepare_cachefile(
    ::Type{D}, input::A, basename::AbstractString
) where {D<:MaxMin, E<:Number, A<:AbstractArray{E, <:Any}}
    mm = MaxMin(input, 10)
    npair = length(mm)
    return (npair, 2, E, mm)
end
