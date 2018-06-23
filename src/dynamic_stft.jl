struct DynCachingStftPsd{
    C<:CachingStftPsd, D<:DynamicWindower
} <: AbstractDynamicSpectrogram{Float64}
    cacher::C
    winput::D
    function DynCachingStftPsd{C, D}(
        cacher::C, winput::D
    ) where {C<:CachingStftPsd, D<:DynamicWindower}
        if ! (basedata(cacher) === basedata(winput))
            throw(ArgumentError("winput and cacher must have same underlying data"))
        end
        new(cacher, winput)
    end
end

function DynCachingStftPsd(
    cacher::C, winput::W
) where {C<:CachingStftPsd, W<:DynamicWindower}
    DynCachingStftPsd{C, W}(cacher, winput)
end

function DynCachingStftPsd(cacher::CachingStftPsd, f_overlap::Real = 0)
    sp = stftpsd(cacher)
    cache_windower = cacher.cacher.winput
    winput = DynamicWindower(
        basedata(cacher),
        sp.stft.fs,
        cacher.base_offset,
        cache_windower.dim,
        f_overlap,
        sp.stft.winput.binsize
    )
    DynCachingStftPsd(cacher, winput)
end

function DynCachingStftPsd(
    input::AbstractVector,
    binsize::Integer,
    fs::Real,
    base_offset::Real = 0,
    f_overlap::Real = 0.8,
    args...;
    kwargs...
)
    DynCachingStftPsd(
        CachingStftPsd(input, binsize, fs, base_offset, args...; kwargs...),
        f_overlap
    )
end

stftpsd(csp::DynCachingStftPsd) = stftpsd(csp.cacher)
frequencies(csp::DynCachingStftPsd) = frequencies(csp.cacher)

function extrema(csp::DynCachingStftPsd)
    freqs = frequencies(stftpsd(csp))
    (freqs[1], freqs[end])
end

fs(a::DynCachingStftPsd) = fs(a.winput)
baselength(a::DynCachingStftPsd) = baselength(a.winput)
start_time(a::DynCachingStftPsd) = start_time(a.winput)

function downsamp_req(
    csp::DynCachingStftPsd,
    xb::Real,
    xe::Real,
    reqpts::Integer,
    exact::Bool = true,
    ::Type{T} = Int
) where T<:Integer
    binsize, overlap, i_begin, i_end, nbase, decno, _ =
        cacher_downsamp_info(csp.cacher.cacher, xb, xe, reqpts)
    if decno > 0
        # defer to CacheAccessor
        (times, ys, was_downsampled) = downsamp_req(
            csp.cacher.cacher, xb, xe, reqpts, exact, T
        )
        freqs = frequencies(csp)
        fwidth = csp.cacher.fwidth
    else
        # make new stft
        binsize, overlap, i_begin, i_end = downsamp_req_window_info(
            csp.winput, xb, xe, reqpts)
        times, wa, was_downsampled = make_windowedarray(
            csp.winput, binsize, overlap, i_begin, i_end
        )
        ft  = stftpsd(csp).stft
        if binsize <= ft.nfft
            downsampler = StftPsd(wa, ft.plan, ft.fs; winfun = ft.winfun)
            st = downsampler.stft
            out  = make_out(st)
            out[:] = 0
            v = st.winput[1]
            GLTimeseries.window_index!(st, 1)
            getindex!(out, st, 1)
            freqs = frequencies(csp)
            fwidth = csp.cacher.fwidth
        else
            downsampler = StftPsd(wa, ft.fs; winfun = ft.winfun)
            freqs = downsampler.stft.frequencies
            fwidth = length(freqs) > 1 ? freqs[2] - freqs[1] : ft.fs / 2
        end
        ys = collect(downsampler)
    end
    y_mat = cat(2, ys...)
    if length(times) > 1
        twidth = times[2] - times[1]
    else
        twidth = binsize / fs(csp)
    end
    (times, (frequencies(csp), y_mat, twidth, fwidth), was_downsampled)
end
