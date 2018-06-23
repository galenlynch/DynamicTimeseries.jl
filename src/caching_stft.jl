struct CachingStftPsd{C<:CacheAccessor} <: AbstractDynamicSpectrogram{Float64}
    cacher::C
    fwidth::Float64
    base_offset::Float64
    function CachingStftPsd{C}(ca::C) where
        {
            C<:CacheAccessor{
                <:DynamicWindower{<:Any, <:Any, <:Any, <:StftPsd},
                <:Any, <:Any, <:Any, <:Any
            }
        }

        sp = basedata(ca.winput)
        base_offset = ca.winput.offset  - (sp.stft.winput.binsize / sp.stft.fs)
        freqs = frequencies(sp)
        fwidth = length(freqs) > 1 ? freqs[2] - freqs[1] : sp.fs / 2
        new(ca, fwidth, base_offset)
    end
end

CachingStftPsd(ca::C) where {C<:CacheAccessor} = CachingStftPsd{C}(ca)

function CachingStftPsd(
    sp::StftPsd, base_offset::Real = 0, args...; kwargs...
)
    new_fs = sp.stft.fs / sp.stft.winput.binsize
    new_offset = base_offset + (sp.stft.winput.binsize / sp.stft.fs)
    CachingStftPsd(
        CacheAccessor(Averager, sp, new_fs, new_offset, args...; kwargs...)
    )
end

function CachingStftPsd(
    input::AbstractVector{<:Number}, binsize::Integer, fs::Real, args...;
    winfun::Function = blackman, kwargs...
)
    CachingStftPsd(
        StftPsd(input, binsize, fs; winfun = winfun), args...; kwargs...
    )
end

stftpsd(csp::CachingStftPsd) = basedata(csp.cacher.winput)
frequencies(csp::CachingStftPsd) = frequencies(stftpsd(csp))
basedata(csp::CachingStftPsd) = basedata(stftpsd(csp))

function extrema(csp::CachingStftPsd)
    freqs = frequencies(csp)
    (freqs[1], freqs[end])
end

fs(csp::CachingStftPsd)  = fs(csp.cacher)
baselength(csp::CachingStftPsd) = baselength(csp.cacher)
start_time(csp::CachingStftPsd) = start_time(csp.cacher)

function downsamp_req(csp::CachingStftPsd, args...)
    (xs, ys, wd) = downsamp_req(csp.cacher, args...)
    y_mat = cat(2, ys...)
    sp = basedata(csp.cacher.winput)
    twidth = length(xs) > 1 ? xs[2] - xs[1] : sp.stft.winput.binsize / sp.stft.fs
    (xs, (frequencies(csp), y_mat, twidth, csp.fwidth), wd)
end
