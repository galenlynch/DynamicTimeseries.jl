abstract type AverageDownsampler{E} <: DynamicDownsampler{E} end

struct AveragingTs{E,W<:DynamicWindower} <: AverageDownsampler{E}
    winput::W
end

function AveragingTs(w::W) where {T<:Number,N,W<:DynamicWindower{<:Any,T,N,<:Any}}
    S = div_type(T)
    E = Array{S, N}
    AveragingTs{E,W}(w)
end

function AveragingTs(w::W) where {T<:Number,W<:DynamicWindower{<:Any,T,1,<:Any}}
    S = div_type(T)
    AveragingTs{S,W}(w)
end

function AveragingTs(
    input::AbstractArray,
    fs::Real,
    offset::Real = 0,
    dim::Integer = 1,
    f_overlap::Real = 0.0,
)
    AveragingTs(DynamicWindower(input, fs, offset, dim, f_overlap))
end

function downsamp_req(ds::AveragingTs{E,<:Any}, xb, xe, npt) where E
    (xs, wa, wd) = downsamp_req(ds.winput, xb, xe, npt)
    avg = Averager(wa)
    ys = collect(avg)
    return xs, ys, wd
end

fs(d::AveragingTs) = fs(d.winput)
baselength(d::AveragingTs) = baselength(d.winput)
start_time(d::AveragingTs) = start_time(d.winput)
