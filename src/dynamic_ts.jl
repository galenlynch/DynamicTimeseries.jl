struct DynamicTs{S<:Number, A<:AbstractVector{S}} <: ExtremaDownsampler{S}
    input::A
    fs::Float64
    offset::Float64
    function DynamicTs{S, A}(
        input::A,
        fs::Float64,
        offset::Float64
    ) where {S<:Number, A<:AbstractVector{S}}
        return new(input, fs, offset)
    end
end
function DynamicTs(
    input::A,
    fs::Real,
    offset::Real = 0
) where {S<:Number, A <: AbstractVector{S}}
    return DynamicTs{S,A}(input, convert(Float64, fs), convert(Float64, offset))
end

function downsamp_req(
    dts::DynamicTs{S,A}, x_start, x_end, reqpoints::Integer, args...
) where {S, A}
    (in_range, i_begin, i_end) = downsamp_range_check(dts, x_start, x_end)
    npt = n_ndx(i_begin, i_end)

    if ! in_range || reqpoints == 0 || npt == 0
        return (Vector{Float64}(), Vector{NTuple{2,S}}(), false)
    end

    # Calculate binsize
    binsize = max(fld(npt, reqpoints), 1)

    was_downsampled = binsize > 1

    # Clip data
    subview = view(dts.input, i_begin:i_end)

    # Create downsampled data
    mm = MaxMin(subview, binsize)
    bin_start_x = ndx_to_t(i_begin, dts.fs, dts.offset)
    xs = ndx_to_t(bin_center(bin_bounds(mm)), dts.fs, bin_start_x)
    ys = collect(mm)
    return (xs::Vector{Float64}, ys::Vector{NTuple{2,S}}, was_downsampled::Bool)
end

baselength(a::DynamicTs) = length(a.input)
start_time(a::DynamicTs) = a.offset
fs(a::DynamicTs) = a.fs
