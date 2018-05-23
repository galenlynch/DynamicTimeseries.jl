"Must implement downsamp_req and duration and extrema"
abstract type DynamicDownsampler{E} end

function extrema(d::DynamicDownsampler)
    (_, ys, _) = downsamp_req(d, duration(d)..., one(Int32), false)
    extrema_red(ys)
end

function downsamp_req(
    ds::DynamicDownsampler, x_start, x_end, reqpoints::AbstractFloat, args...
)
    return downsamp_req(ds, x_start, x_end, floor(Int, reqpoints), args...)
end

function downsamp_range_check(
    dts::DynamicDownsampler, x_start, x_end, ::Type{T} = Int32
) where {T<:Integer}
    x_start <= x_end || throw(ArgumentError("x_start must be before x_end"))
    nin = T(length(dts.input))
    # Find bounding indices
    i_begin = t_to_ndx(x_start, dts.fs, dts.offset, T)
    i_end = t_to_ndx(x_end, dts.fs, dts.offset, T)
    in_range = i_begin <= nin && i_end >= 1
    i_begin_clipped = clip_ndx(i_begin, nin)
    i_end_clipped = clip_ndx(i_end, nin)
    return (in_range, i_begin_clipped, i_end_clipped)
end
