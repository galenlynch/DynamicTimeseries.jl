struct DynamicSpectrogram{A<:AbstractVector} <: DynamicDownsampler{Tuple{Vector{Float64}, Array{Float64, 2}}}
    input::A
    fs::Float64
    offset::Float64
    window::Vector{Float64}
    overlap::Float64
    function DynamicSpectrogram{A}(
        input::A,
        fs::Float64,
        offset::Float64,
        window::Vector{Float64},
        overlap::Float64
    ) where {A <: AbstractVector}
        @assert fs > 0 "fs must be greater than zero"
        @assert 0 <= overlap < 1 "overlap must be in interval [0, 1)"
        return new(input, fs, offset, window, overlap)
    end
end
function DynamicSpectrogram(
    input::A,
    fs::Real,
    offset::Real = 0,
    window::Vector{Float64} = hanning(512),
    overlap::Float64 = 0.8
) where {A<:AbstractVector}
    return DynamicSpectrogram{A}(
        input,
        convert(Float64, fs),
        convert(Float64, offset),
        window,
        convert(Float64, overlap)
    )
end

duration(d::DynamicSpectrogram) = duration(length(d.input), d.fs, d.offset)

function extrema(d::DynamicSpectrogram)
    freqs = rfftfreq(length(d.window), d.fs)
    return (freqs[1], freqs[end])
end

"""
    downsamp_req(ds, xb, xe, npt; windowfun)

Calculates the spectrogram for data in the range between xb and xe. Expands the
selected range so the center of the first time bin is around xb, and the
center of the last time bin is around xe.
"""
function downsamp_req(
    ds::DynamicSpectrogram, xb, xe, npt::Integer, args...;
    windowfun::Function = hanning
)
    # Figure out what was asked of us
    nin = length(ds.input)
    (ib, ie) = clip_ndx.(t_to_ndx.([xb, xe], ds.fs, ds.offset), nin)
    win_l = length(ds.window)
    nsel = n_ndx(ib, ie)
    noverlap = floor(Int, ds.overlap * win_l)

    # let:
    # N := number of points in the selected region
    # l := number of points in a window
    # w := number of windows
    # d := number of overlapping samples
    # Then we have N = (w - 1) * l - (w - 1) * d
    # => w = (N + l - d) / (l - d)
    win_nbase = div(nin + win_l - noverlap, win_l - noverlap)

    # Check if we have too many windows
    if win_nbase > npt
        # Too many windows: increase their size

        # f := fraction of window to overlap (=> d = fl)
        # Substitute d = fl into the above, yielding:
        # l = N / (w - f * (w - 1) - 1)
        win_l_target = cld(nsel, npt - floor(Int, ds.overlap * (npt - 1)) - 1)
        was_downsamped = true
    else
        win_l_target = win_l
        was_downsamped = false
    end

    # Check if we have enough data for our windows
    if nsel < win_l_target
        win_l_final = nsel
        npt_overlap_final = 0
    else
        win_l_final = win_l_target
        npt_overlap_final = floor(Int, ds.overlap * win_l_final)
    end

    # Check if our window function is the right size
    if length(ds.window) != win_l_final
        window = windowfun # Pass the function itself
    else
        window = ds.window
    end

    # Expand selected data for the bins on the edge
    half_w = div(win_l_final, 2)
    ib_ex = max(1, ib - half_w)
    ie_ex = min(nin, ie + half_w)

    # Select the portion of the signal in range
    sel_sig = ds.input[ib_ex:ie_ex]

    S = spectrogram(
        sel_sig, win_l_final, npt_overlap_final;
        fs = ds.fs, window = window
    )

    first_x = ndx_to_t(ib_ex, ds.fs, ds.offset)
    times = S.time + first_x

    t_width = length(times) > 1 ? times[2] - times[1] : win_l_final / ds.fs
    f_width = length(S.freq) > 1 ? S.freq[2] - S.freq[1] : ds.fs / 2

    return (times, (collect(S.freq), S.power, t_width, f_width), was_downsamped)
end
