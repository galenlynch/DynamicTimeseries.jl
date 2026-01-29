struct DynCachingStftPsd{C<:CachingStftPsd,W<:DynamicWindower} <:
       AbstractDynamicSpectrogram{Float64}
    cachepsd::C
    winput::W
    function DynCachingStftPsd{C,W}(
        cachepsd::C,
        winput::W,
    ) where {C<:CachingStftPsd,W<:DynamicWindower}
        if ! (basedata(cachepsd) === basedata(winput))
            throw(ArgumentError("winput and cachepsd must have same underlying data"))
        end
        new(cachepsd, winput)
    end
end

function DynCachingStftPsd(
    cachepsd::C,
    winput::W,
) where {C<:CachingStftPsd,W<:DynamicWindower}
    DynCachingStftPsd{C,W}(cachepsd, winput)
end

function DynCachingStftPsd(cachepsd::CachingStftPsd, f_overlap::Real = 0)
    sp = stftpsd(cachepsd)
    cache_windower = cache_accessor(cachepsd).winput
    winput = DynamicWindower(
        basedata(cachepsd),
        sp.stft.fs,
        cachepsd.base_offset,
        f_overlap,
        sp.stft.winput.binsize,
    )
    DynCachingStftPsd(cachepsd, winput)
end

function DynCachingStftPsd(
    input::AbstractVector,
    binsize::Integer,
    fs::Real,
    winfun::Union{Function,Missing} = blackman,
    base_offset::Real = 0,
    f_overlap::Real = 0.8,
    args...;
    kwargs...,
)
    winfun = ismissing(winfun) ? blackman : winfun
    DynCachingStftPsd(
        CachingStftPsd(input, binsize, fs, 0, winfun, base_offset, args...; kwargs...),
        f_overlap,
    )
end

stftpsd(csp::DynCachingStftPsd) = stftpsd(csp.cachepsd)
frequencies(csp::DynCachingStftPsd) = frequencies(csp.cachepsd)

function extrema(csp::DynCachingStftPsd)
    freqs = frequencies(stftpsd(csp))
    (freqs[1], freqs[end])
end

fs(a::DynCachingStftPsd) = fs(a.winput)
baselength(a::DynCachingStftPsd) = baselength(a.winput)
start_time(a::DynCachingStftPsd) = start_time(a.winput)
basedata(a::DynCachingStftPsd) = basedata(a.winput)

function cache_accessor(a::DynCachingStftPsd{C,W}) where {W,E,C<:CachingStftPsd{E}}
    cache_accessor(a.cachepsd)::E
end

function downsamp_req(
    csp::DynCachingStftPsd{C},
    xb::Real,
    xe::Real,
    reqpts::Integer,
    exact::Bool = false,
    ::Type{T} = Int,
) where {T<:Integer,C}
    base_cacher = cache_accessor(csp)
    base_psd = csp.cachepsd::C
    binsize, overlap, i_begin, i_end, nbase, decno, _ =
        cacher_downsamp_info(base_cacher, xb, xe, reqpts)
    if decno > 0
        # defer to CacheAccessor
        (times, ys, was_downsampled) = downsamp_req(base_cacher, xb, xe, reqpts, exact, T)
        freqs = frequencies(csp)
        fwidth = base_psd.fwidth
        y_len = isempty(ys) ? 0 : length(ys[1])
        ny = length(ys)
        @compat y_mat = Array{eltype(eltype(ys)),2}(undef, y_len, ny)
        @inbounds @simd for i = 1:ny
            y_mat[:, i] = ys[i]
        end
    else
        # make new stft
        binsize, overlap, i_begin, i_end =
            downsamp_req_window_info(csp.winput, xb, xe, reqpts)
        times, wa, was_downsampled =
            make_windowedarray(csp.winput, binsize, overlap, i_begin, i_end)
        ft = stftpsd(csp).stft
        if binsize <= ft.nfft
            downsampler = StftPsd(wa, ft.plan, ft.fs, ft.winfun)
            st = downsampler.stft
            out = make_out(st)
            out[:] .= 0
            v = st.winput[1]
            GLTimeseries.window_index!(st, 1)
            getindex!(out, st, 1)
            freqs = frequencies(csp)
            fwidth = base_psd.fwidth
        else
            downsampler = StftPsd(wa, ft.fs, ft.winfun)
            freqs = downsampler.stft.frequencies
            fwidth = length(freqs) > 1 ? freqs[2] - freqs[1] : ft.fs / 2
        end
        @compat y_mat =
            Array{Float64,2}(undef, prod(elsize(downsampler)), length(downsampler))
        copy!(y_mat, downsampler)
    end
    if length(times) > 1
        twidth = times[2] - times[1]
    else
        twidth = binsize / fs(csp)
    end
    (times, (freqs, y_mat, twidth, fwidth), was_downsampled)
end
