struct CachingStftPsd{C<:CacheAccessor}
    cacher::C
end

function CachingStftPsd(sp::StftPsd, args...; kwargs...)
    new_fs = sp.stft.fs / sp.stft.winput.binsize
    CachingStftPsd(
        CacheAccessor(Averager, sp, new_fs, args...; kwargs...)
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

function extrema(csp::CachingStftPsd)
    freqs = frequencies(basedata(csp.cacher.winput))
    (freqs[1], freqs[end])
end

fs(csp::CachingStftPsd)  = fs(csp.cacher)
downsamp_req(csp::CachingStftPsd, args...) = downsamp_req(csp.cacher, args...)
baselength(csp::CachingStftPsd) = baselength(csp.cacher)
start_time(csp::CachingStftPsd) = start_time(csp.cacher)
