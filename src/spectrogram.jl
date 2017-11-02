struct DynamicSpectrogram{A<:AbstractVector} <: DynamicDownsampler
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
    window::Vector{Float64} = hanning(1024),
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

function downsamp_req(ds::DynamicSpectrogram, xb, xe, npt::Integer; windowfun::Function = hanning)
    # test
    nin = length(ds.input)
    (ib, ie) = clip_ndx.(t_to_ndx.([xb, xe], ds.fs, ds.offset), nin)
    nsel = n_ndx(ib, ie)
    win_l = length(ds.window)
    win_nbase = cld(nsel, ceil(Int, (1 - ds.overlap) * win_l)) # Number of windows that would be used

    # Check if we have too many windows
    if win_nbase > npt
        # Too many windows: increase their size
        win_l_target = cld(nsel, ceil(Int, (1 - ds.overlap) * npt))
    else
        win_l_target = win_l
    end

    # Check if we have enough data for our windows
    if nsel < win_l_target
        win_l_final = nsel
        npt_overlap = 0
    else
        win_l_final = win_l_target
        npt_overlap = floor(Int, ds.overlap * win_l_final)
    end

    # Check if our window function is the right size
    if length(ds.window) != win_l_final
        window = hanning(win_l_final)
    else
        window = ds.window
    end

    # Select the portion of the signal in range
    sel_sig = ds.input[ib:ie]
    S = spectrogram(sel_sig, win_l_final, npt_overlap; fs = ds.fs, window = window)
    first_x = ndx_to_t(ib, ds.fs, ds.offset)
    times = S.time + first_x
    return (times, (collect(S.freq), S.power))
end
