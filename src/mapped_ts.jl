struct MappedDynamicDownsampler{E, D<:DynamicDownsampler{E}} <: AbstractDynamicDownsampler{E}
    downsampler::D
    fmap::Function
end
function MappedDynamicDownsampler(d::D, fmap::Function) where {E, D<:DynamicDownsampler{E}}
    return MappedDynamicDownsampler{E,D}(d, fmap)
end

function downsamp_req(
    mdds::D, xb, xe, npts::Integer, exact::Bool = false, ::Type{T} = Int32
) where {D<:MappedDynamicDownsampler, T<:Integer}
    (xs, ys, was_downsampled) = downsamp_req(mdds.downsampler, xb, xe, npts, exact, T)
    mys = mdds.fmap(ys)
    return (xs, mys, was_downsampled)
end

duration(d::MappedDynamicDownsampler) = duration(d.downsampler)
fs(d::MappedDynamicDownsampler) = fs(d.downsampler)
baselength(d::MappedDynamicDownsampler) = baselength(d.downsampler)
start_time(d::MappedDynamicDownsampler) = start_time(d.downsampler)

make_shifter(shift) = (x) -> shift_extrema!(shift, x)

function shift_extrema(shift, ys::T) where {S, T<:NTuple{2,S}}
    (ys[1] + shift, ys[2] + shift)
end
shift_extrema(shift, y::Number) = y + shift

function shift_extrema!(shift, dest::B, ys::A) where
    {S, A<:AbstractVector{S}, B<:AbstractVector{S}}
    for (i, tup) in enumerate(ys)
        dest[i] = shift_extrema(shift, ys[i])
    end
end

function shift_extrema!(shift, ys)
    shift_extrema!(shift, ys, ys)
    return ys
end

function shift_extrema(shift, ys::AbstractVector)
    shifted = similar(ys)
    shift_extrema!(shift, shifted, ys)
    return shifted
end
