struct DynamicSpectrogram{A<:AbstractVector} <: DynamicDownsampler
    input::A
    fs::Float64
    offset::Float64
    window::Vector{Float64}
end
function DynamicSpectrogram(
    input::A,
    fs::Real,
    offset::Real = 0,
    window::Vector{Float64} = hanning(512)
) where {A<:AbstractVector}
    return DynamicSpectrogram{A}(input, convert(Float64, fs), convert(Float64, offset), window)
end

function downsamp_req(ds::DynamicSpectrogram, xb::Real, xe::Real, npt::Integer)
    nin = length(ds.input)
    (ib, ie) = clipind.(x_to_ndx.((xb, xe), ds.fs, ds.offset), nin)
    nsel = n_ndx(ib, ie)
    win_l = length(ds.window)
    if nsel >= win_l
        win_nbase = cld(nsel, win_l)
        if win_nbase > npt
            # Too many windows: increase their size
            win_l_final = 2 * cld(nsel, npt)
            window = hanning(win_l_final)
            overlap_f = 0.5
        else
            # Too few windows: increase their overlap_f
            window = ds.window
            win_l_final = win_l
            npt_per_win = cld(npt, win_nbase)
            overlap_f = min(0.999, 1 - (1 / npt_per_win))
        end
        npt_overlap = convert(Int, floor(overlap_f * win_l_final))
    else
        window = ds.window
        win_l_final = win_l
        npt_overlap = 0
    end
    S = spectrogram(ds.input, win_l_final, npt_overlap; fs = ds.fs, window = window)
    first_x = ndx_to_x(ib, ds.fs, ds.offset)
    times = S.time + first_x
    return (S.power, collect(S.freq), times)
end
function downsamp_req(ds::DynamicSpectrogram, x_start::Real, x_end::Real, reqpoints::AbstractFloat)
    return downsamp_req(ds, x_start, x_end, convert(Int, floor(reqpoints)))
end
