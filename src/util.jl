function extent(series::A) where {T<:Number, A<:AbstractVector{T}}
    extremum = extrema(series)
    return extremum[2] - extremum[1]
end
function extent(series::A) where {T<:AbstractVector, A<:AbstractVector{T}}
    return extent.(series)
end
