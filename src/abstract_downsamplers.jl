"""
Must implement
array interface
binsize
bin_bounds
downsamp_reduce
downsamp_reduce_cache
"""
abstract type Downsampler{T, N} <: AbstractArray{T, N} end

function downsamp_reduce_cache(::Type{D}, args...) where D<:Downsampler
    downsamp_reduce(D, args...)
end

"""
Must implement:
downsamp_req
extrema
fs
baselength
start_time
"""
abstract type AbstractDynamicDownsampler{E} end

function downsamp_req(
    ds::AbstractDynamicDownsampler, x_start, x_end, reqpoints::AbstractFloat, args...
)
    return downsamp_req(ds, x_start, x_end, floor(Int, reqpoints), args...)
end

function downsamp_range_check(
    dts::AbstractDynamicDownsampler, x_start::Real, x_end::Real, ::Type{T} = Int32
) where {T<:Integer}
    x_start <= x_end || throw(ArgumentError("x_start must be before x_end"))
    nin = T(baselength(dts))
    # Find bounding indices
    d_fs = fs(dts)
    d_ss = start_time(dts)
    i_begin = t_to_ndx(x_start, d_fs, d_ss, T)
    i_end = t_to_last_ndx(x_end, d_fs, d_ss, T)
    in_range = i_begin <= nin && i_end >= 1
    i_begin_clipped = clip_ndx(i_begin, nin)
    i_end_clipped = clip_ndx(i_end, nin)
    return (in_range, i_begin_clipped::T, i_end_clipped::T)
end

duration(a::AbstractDynamicDownsampler) = baselength(a) / fs(a)
time_interval(a::AbstractDynamicDownsampler) = (start_time(a), start_time(a) + duration(a))
