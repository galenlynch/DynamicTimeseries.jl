struct DynamicTs{S<:Number, W<:DynamicWindower{<:Any,S,1,<:Any}} <: ExtremaDownsampler{S}
    winput::W
end

function DynamicTs(a::W) where {T<:Number,W<:DynamicWindower{<:Any,T,1,<:Any}}
    DynamicTs{T,W}(a)
end

function DynamicTs(
    input::A,
    fs::Real,
    offset::Real = 0
) where {S<:Number, A <: AbstractVector{S}}
    winput = DynamicWindower(input, fs, offset)
    DynamicTs(winput)
end

function downsamp_req(
    dts::DynamicTs{S,A}, x_start, x_end, reqpoints::Integer, args...
) where {S, A}
    (xs, wa, wd) = downsamp_req(dts.winput, x_start, x_end, reqpoints)
    mm = MaxMin(wa)
    ys = collect(mm)
    return (xs::Vector{Float64}, ys::Vector{NTuple{2,S}}, wd)
end

baselength(a::DynamicTs) = baselength(a.winput)
start_time(a::DynamicTs) = start_time(a.winput)
fs(a::DynamicTs) = fs(a.winput)
