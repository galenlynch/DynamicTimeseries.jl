struct Stft{
    E, W<:WindowedArray, P<:Base.DFT.Plan
} <: Downsampler{E, 1}
    winput::W
    fs::Float64
    win::Vector{Float64}
    winfun::Function # must be Integer -> Vector{Float64}
    r::Float64
    r_temp::Vector{Float64}
    plan::P
    winbuf::Vector{Float64}
    frequencies::Vector{Float64}
    nfft::Int
    nout::Int
    function Stft{E, W, P}(
        winput::W, plan::P, fs::Float64, winfun::Function
    ) where {E, W<:WindowedArray, P<:Base.DFT.Plan}

        nfft = length(plan)
        nout = div(nfft, 2) + 1

        win = winfun(winput.binsize)::Vector{Float64}

        norm2 = sum(abs2, win)
        r = fs * norm2
        r_temp = Vector{Float64}(1)

        fftbuf = Vector{Complex{Float64}}(nout)
        winbuf = zeros(Float64, nfft)

        frequencies = collect(rfftfreq(nfft, fs))

        new(
            winput,
            fs,
            win,
            winfun,
            r,
            r_temp,
            plan,
            winbuf,
            frequencies,
            nfft,
            nout
        )
    end
end

function Stft(
    winput::W, plan::P, fs::Real; winfun::Function = blackman
) where {W<:WindowedArray{<:Any, <:Any, 1, <:Any}, P<:Base.DFT.Plan}
    Stft{Vector{Complex{Float64}}, W, P}(
        winput, plan, convert(Float64, fs), winfun
    )
end

function Stft(
    winput::W, fs::Real = 1, args...; kwargs...
) where W<:WindowedArray{<:Any, <:Any, 1, <:Any}
    nfft = nextfastfft(winput.binsize)
    plan = plan_rfft(Vector{Float64}(nfft))
    Stft(winput, plan, fs, args...; kwargs...)
end

function Stft(
    input::AbstractVector,
    binsize::Integer,
    fs::Real,
    args...;
    overlap::Integer = 0,
    kwargs...
)
    winput = WindowedArray(input, binsize, 1, overlap)
    Stft(winput, fs, args...; kwargs...)
end

frequencies(a::Stft) = a.frequencies

function getindex(a::Stft, i::Integer)
    out = make_out(a)
    getindex!(out, a, i)
end

function getindex!(out, a::Stft, i::Integer)
    @boundscheck checkbounds(a, i)
    window_index!(a, i)
    A_mul_B!(out, a.plan, a.winbuf)
end

make_out(a::Stft) = Vector{Complex{Float64}}(a.nout)

function window_index!(a::Stft, i::Integer)
    v = a.winput[i]
    nv = length(v)
    full_bin = length(a.win) == nv

    if full_bin
        win = a.win
        a.r_temp[1] = a.r
    else
        win = a.winfun(nv)::Vector{Float64}
        a.r_temp[1] = a.fs * sum(abs2, win)
    end
    @simd for i in 1:nv
        @inbounds a.winbuf[i] = win[i] * v[i]
    end
    @simd for i in (nv + 1):a.winput.binsize
        @inbounds a.winbuf[i] = 0
    end
end

function concrete_type(::Type{<:Stft}, ::Type{W}) where {W<:WindowedArray}
    return Stft{
        Vector{Float64},
        W,
        rFFTWPlan{Float64, Base.DFT.FFTW.FORWARD, false, 1}
    }
end
length(a::Stft) = length(a.winput)
size(a::Stft) = size(a.winput)

setindex!(::Stft) = throw(ReadOnlyMemoryError())

Base.IndexStyle(::Type{Stft{<:Any, W, <:Any}}) where W = IndexStyle(W)

binsize(a::Stft) = binsize(a.winput)

bin_bounds(a::Stft, args...) = bin_bounds(a.winput, args...)
bin_bounds(i::Integer, a::Stft) = bin_bounds(i, a.winput)

struct StftPsd{E, S<:Stft} <: Downsampler{E, 1}
    stft::Stft
    fftbuf::Vector{Complex{Float64}}
    function StftPsd{E, S}(stft::Stft) where {E, S<:Stft}
        fftbuf = make_out(stft)
        new(stft, fftbuf)
    end
end
StftPsd(stft::S) where S<:Stft = StftPsd{Vector{Float64}, S}(stft)
StftPsd(args...; kwargs...) = StftPsd(Stft(args...; kwargs...))

make_out(a::StftPsd) = Vector{Float64}(a.stft.nout)

function getindex(s::StftPsd, i::Integer)
    dest = make_out(s)
    getindex!(dest, s, i)
end

function getindex!(out, s::StftPsd, i::Integer)
    @boundscheck checkbounds(s, i)
    getindex!(s.fftbuf, s.stft, i)
    fft2pow!(out, s.fftbuf, s.stft.nfft, s.stft.r_temp[1])
end

frequencies(s::StftPsd) = frequencies(s.stft)

function fft2pow!(
    dest::Vector{T}, s_fft::Vector{Complex{T}}, nfft::Integer, r::Real
) where T
    n = length(s_fft)
    n <= length(dest) || throw(ArgumentError("dest must be at least as long as s_fft"))
    m1 = convert(T, 1/r)
    m2 = 2 * m1
    @inbounds dest[1] = m1 * abs2(s_fft[1])
    @simd for i = 2:(n - 1)
        @inbounds dest[i] = m2 * abs2(s_fft[i])
    end
    @inbounds dest[n] = ifelse(iseven(nfft), m1, m2) * abs2(s_fft[n])
    dest
end

function concrete_type(::Type{<:StftPsd}, ::Type{W}) where {W<:WindowedArray}
    return StftPsd{Vector{Float64}, concrete_type(Stft, W)}
end
length(a::StftPsd) = length(a.stft)
size(a::StftPsd) = size(a.stft)

setindex!(::StftPsd) = throw(ReadOnlyMemoryError())

Base.IndexStyle(::Type{StftPsd{<:Any, S}}) where S<:Stft = IndexStyle(S)

binsize(a::StftPsd) = binsize(a.winput)

bin_bounds(a::StftPsd, args...) = bin_bounds(a.stft, args...)
bin_bounds(i::Integer, a::StftPsd) = bin_bounds(i, a.stft)
