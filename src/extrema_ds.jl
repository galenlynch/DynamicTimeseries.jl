abstract type ExtremaDownsampler{E} <: DynamicDownsampler{NTuple{2, E}} end

function extrema(d::ExtremaDownsampler)
    (_, ys, _) = downsamp_req(d, time_interval(d)..., one(Int32), false)
    extrema_red(ys)
end
