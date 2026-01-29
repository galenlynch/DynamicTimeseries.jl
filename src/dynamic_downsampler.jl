struct DynamicDownsampler{E,D<:Downsampler,W<:DynamicWindower} <:
       AbstractDynamicDownsampler{E}
    downsampler::Type{D}
    winput::W
end

function DynamicDownsampler(
    ::Type{D},
    winput::W,
) where {D<:Downsampler,A,W<:DynamicWindower{<:Any,<:Any,<:Any,A}}
    DynamicDownsampler{eltype_preview(D, A),D,W}(D, winput)
end

function DynamicDownsampler(::Type{D}, input::AbstractArray, args...) where {D<:Downsampler}
    DynamicDownsampler(D, DynamicWindower(input, args...))
end

function downsamp_req(
    dd::D,
    xb,
    xe,
    npt,
    args...,
) where {E,C,D<:DynamicDownsampler{E,C,<:Any}}
    (xs, wa, wd) = downsamp_req(dd.winput, xb, xe, npt)
    ds = C(wa)
    ys = collect(ds)
    return xs, ys::Vector{E}, wd
end

baselength(a::DynamicDownsampler) = baslength(a.winput)
start_time(a::DynamicDownsampler) = start_time(a.winput)
fs(a::DynamicDownsampler) = fs(a.winput)
