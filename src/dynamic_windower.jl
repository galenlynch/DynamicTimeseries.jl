struct DynamicWindower{E,T,N,A<:AbstractArray{T,N}} <: AbstractDynamicDownsampler{E}
    input::A
    fs::Float64
    offset::Float64
    f_overlap::Float64
    wmin::Int
    function DynamicWindower{E,T,N,A}(
        input::A,
        fs::Float64,
        offset::Float64 = 0.0,
        f_overlap::Float64 = 0.0,
        wmin::Int = 1,
    ) where {E,T,N,A<:AbstractArray{T,N}}
        validate_dynamic_windower_args(N, f_overlap, wmin, fs)
        return new(input, fs, offset, f_overlap, wmin)
    end
end

function DynamicWindower(
    input::A,
    fs::Real,
    offset::Real = 0,
    f_overlap::Real = 0,
    wmin::Integer = 1,
) where {T,N,A<:AbstractArray{T,N}}
    v = view_trailing_slice(input, 1:0)
    return DynamicWindower{typeof(v),T,N,A}(
        input,
        convert(Float64, fs),
        convert(Float64, offset),
        convert(Float64, f_overlap),
        convert(Int, wmin),
    )
end

function validate_dynamic_windower_args(
    N::Integer,
    f_overlap::AbstractFloat,
    wmin::Integer,
    fs::Real,
)
    if fs <= 0
        throw(ArgumentError("fs must be greater than zero"))
    end
    if f_overlap < 0 || f_overlap >= 1
        throw(ArgumentError("f_overlap must be in interval [0, 1), but it is $f_overlap"))
    end
    if wmin < 1
        throw(ArgumentError("wmin must be at least one"))
    end
end

function windowtype(::Type{W}) where {E,T,N,A,W<:DynamicWindower{E,T,N,A}}
    WindowedArray{E,T,N,A}
end

fs(d::DynamicWindower) = d.fs
baselength(d::DynamicWindower) = length(d.input)
basedata(d::DynamicWindower) = d.input
start_time(d::DynamicWindower) = d.offset

"""
    downsamp_req(ds, xb, xe, npt; windowfun)

    Calculates the spectrogram for data in the range between xb and xe. Expands the
    selected range so the center of the first time bin is around xb, and the
    center of the last time bin is around xe.
    """
function downsamp_req(ds::DynamicWindower, xb, xe, npt::Integer)
    win_l, overlap, ib_ex, ie_ex = downsamp_req_window_info(ds, xb, xe, npt)
    make_windowedarray(ds, win_l, overlap, ib_ex, ie_ex)
end

#TODO type unstable in a bad way
function make_windowedarray(
    ds::DynamicWindower{<:Any,<:Any,N,<:Any},
    win_l,
    overlap,
    ib_ex,
    ie_ex,
) where {N}
    v = view_trailing_slice(ds.input, ib_ex:ie_ex)
    wa = WindowedArray(v, win_l, overlap)

    first_t = ndx_to_t(ib_ex, fs(ds), ds.offset)
    times = ndx_to_t(bin_center(bin_bounds(wa)), fs(ds), first_t)

    was_downsampled = win_l > 1

    return times, wa, was_downsampled
end

function downsamp_req_window_info(ds::DynamicWindower, xb, xe, npt::S) where {S}
    (in_range, ib, ie) = downsamp_range_check(ds, xb, xe, S)
    nsel = n_ndx(ib, ie)
    if ! in_range || npt <= 0 || nsel == 0
        # Bail out early with empty result
        win_l = one(S)
        ib_ex = one(S)
        ie_ex = zero(S)
        overlap = zero(S)
    else
        nin = S(baselength(ds))
        win_l_target = size_windows_expanded(npt, nsel, nin, ds.f_overlap, S(ds.wmin))::S

        if win_l_target > nsel
            win_l = nsel::S
            npt_adj = one(S)
        else
            win_l = win_l_target
            npt_adj = npt
        end

        overlap = floor(S, ds.f_overlap * win_l)

        if npt_adj > 1
            half_w = div(win_l, S(2))
            # Expand start index so the first bin is roughly centered at xb
            ib_ex = max(one(S), ib - half_w)

            # Now some math. If s is the start index of the whole range, then
            # the start index for the nth bin, sn, is:
            # sn = s + (l - d)(n - 1)
            # where l is the window length, and d is the overlap
            #
            # The center of the nth bin, cn, is therefore
            # cn = sn + (l - 1) / 2
            # get the last center around xe, we solve for n.
            # if ie is the index of the last bin in the range, and w is the
            # number of windows, then from the above:
            # ie = cn = sn + (l - 1) / 2 = s + (l-d)(w-1) + (l-1)/2
            # so w = (ie - s + l - (l - 1) / 2) / (l - d)
            # where w is the number of windows
            npt_out = ceil(
                S,
                (ie - ib_ex + win_l - (win_l - one(S)) / S(2)) / (win_l - overlap),
            )::S

            # Finally calculate the last bound, p
            # p = s + w*l - (w - 1)*d - 1
            ie_ex =
                min(nin, ib_ex + npt_out * win_l - (npt_out - one(S)) * overlap - one(S))::S
        else
            ib_ex = ib
            ie_ex = ie
        end
    end
    return win_l::S, overlap::S, ib_ex::S, ie_ex::S
end

function size_windows_expanded(
    n_win_max::T,
    n_point::T,
    win_l_max::T = typemax(T),
    overlap_frac::Real = 0,
    win_l_min::T = one(T),
) where {T<:Integer}
    n_win_max == 1 && return n_point
    n_overlap = floor(T, overlap_frac * win_l_min)
    if n_overlap == win_l_min
        throw(ArgumentException("Overlap must be smaller than window"))
    end

    # let:
    # N := number of points in the selected region
    # l := number of points in a window
    # d := number of overlapping samples in each window
    # w := number of windows
    # Then we have N = (w - 1) * l - (w - 1) * d
    # => w = (N + l - d) / (l - d)

    # number of windows with minimum window length
    n_win = div(n_point + win_l_min - n_overlap, win_l_min - n_overlap)::T

    if n_win > n_win_max
        # Too many windows: increase their size

        # f := fraction of window to overlap (=> d = fl)
        # Substitute d = fl into the above, yielding:
        # l = N / (w - f * (w - 1) - 1)
        win_l =
            cld(n_point, n_win_max - floor(T, overlap_frac * (n_win_max - one(T))) - one(T))
    else
        win_l = win_l_min
    end
    win_l_final = win_l > win_l_max ? win_l_max : win_l
    return win_l_final::T
end

function size_windows_expanded(
    nwm::Integer,
    np::Integer,
    wlmax::Integer,
    of::Real,
    wlmin::Integer,
)
    (nwm, np, wlmax, wlmin) = promote(nwm, np, wlmax, wlmin)
    size_windows_expanded(nwm, np, wlmax, of, wlmin)
end
