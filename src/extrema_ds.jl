const ExtremaDynDownsampler{E,W} = DynamicDownsampler{NTuple{2,E},MaxMin,W}

function extrema(d::ExtremaDynDownsampler)
    (_, ys, _) = downsamp_req(d, time_interval(d)..., one(Int32), false)
    extrema_red(ys)
end
