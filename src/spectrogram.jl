struct DynamicSpectrogram{
    T<:AbstractFloat, W<:DynamicWindower{<:Any,<:Number,1,<:Any}
} <: AverageDownsampler{Tuple{Vector{Float64}, Array{T, 2}}}
    winput::W
    window::Vector{Float64}
    function DynamicSpectrogram{T,W}(
        winput::W, window::Vector{Float64}
    ) where {T<:AbstractFloat,W<:DynamicWindower{<:Any,<:Number,1,<:Any}}
        if length(window) != winput.wmin
            throw(ArgumentError("wmin must match window length"))
        end
        new(winput, window)
    end
end

function DynamicSpectrogram(
    winput::W, window::Vector{Float64} = hanning(512)
) where {T,W<:DynamicWindower{<:Any,T,1,<:Any}}
    S = div_type(T)
    DynamicSpectrogram{S,W}(winput, window)
end

function DynamicSpectrogram(
    input::A,
    fs::Real,
    offset::Real = 0,
    overlap::Float64 = 0.8,
    window::Vector{Float64} = hanning(512),
) where {T<:Real, A<:AbstractVector{T}}
    DynamicSpectrogram(
        DynamicWindower(input, fs, offset, 1, overlap, length(window)),
        window
    )
end

"""
    downsamp_req(ds, xb, xe, npt; windowfun)

Calculates the spectrogram for data in the range between xb and xe. Expands the
selected range so the center of the first time bin is around xb, and the
center of the last time bin is around xe.
"""
function downsamp_req(
    ds::DynamicSpectrogram{T, <:Any}, xb, xe, npt::Integer, args...;
    windowfun::Function = hanning
) where {T<:AbstractFloat}
    # Figure out what was asked of us
    win_l, overlap, ib_ex, ie_ex = downsamp_req_window_info(ds.winput, xb, xe, npt)


    v = view(basedata(ds.winput), ib_ex:ie_ex)
    # Select the portion of the signal in range

    window = win_l == length(ds.window) ? ds.window : hanning

    srate = fs(ds.winput)

    S = spectrogram(
        v, win_l, overlap;
        fs = srate, window = window
    )

    first_x = ndx_to_t(ib_ex, fs(ds.winput), start_time(ds.winput))
    times = S.time + first_x

    t_width = length(times) > 1 ? times[2] - times[1] : win_l / srate
    f_width = length(S.freq) > 1 ? S.freq[2] - S.freq[1] : srate / 2

    return (
        times::StepRangeLen{
            Float64,
            Base.TwicePrecision{Float64},
            Base.TwicePrecision{Float64}
        },
        (
            collect(S.freq)::Vector{Float64},
            S.power::Array{T, 2},
            t_width::Float64,
            f_width::Float64
        ),
        false
    )
end

fs(d::DynamicSpectrogram) = fs(d.winput)
baselength(d::DynamicSpectrogram) = baselength(d.winput)
start_time(d::DynamicSpectrogram) = start_time(d.winput)

function extrema(d::DynamicSpectrogram)
    freqs = rfftfreq(length(d.window), fs(d.winput))
    return (freqs[1], freqs[end])
end
