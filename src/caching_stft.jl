struct CachingStftPsd{C<:CacheAccessor} <: AbstractDynamicSpectrogram{Float64}
    cacher::C
    freqs::Vector{Float64}
    fwidth::Float64
    function CachingStftPsd{C}(ca::C) where
        {
            C<:CacheAccessor{
                <:DynamicWindower{<:Any, <:Any, <:Any, <:StftPsd},
                <:Any, <:Any, <:Any, <:Any
            }
        }
        sp = basedata(ca.winput)
        freqs = frequencies(sp)
        fwidth = length(freqs) > 1 ? freqs[2] - freqs[1] : sp.fs / 2
        new(ca, freqs, fwidth)
    end
end

CachingStftPsd(ca::C) where {C<:CacheAccessor} = CachingStftPsd{C}(ca)

function CachingStftPsd(sp::StftPsd, args...; kwargs...)
    new_fs = sp.stft.fs / sp.stft.winput.binsize
    CachingStftPsd(CacheAccessor(Averager, sp, new_fs, args...; kwargs...))
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
baselength(csp::CachingStftPsd) = baselength(csp.cacher)
start_time(csp::CachingStftPsd) = start_time(csp.cacher)

function downsamp_req(csp::CachingStftPsd, args...)
    (xs, ys, wd) = downsamp_req(csp.cacher, args...)
    y_mat = cat(2, ys...)
    sp = basedata(csp.cacher.winput)
    twidth = length(xs) > 1 ? xs[2] - xs[1] : sp.stft.winput.binsize / sp.stft.fs
    (xs, (csp.freqs, y_mat, twidth, csp.fwidth), wd)
end
