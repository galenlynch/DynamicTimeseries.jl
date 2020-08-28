struct Stft{
    E, W<:WindowedArray, P<:Plan
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
    ) where {E, W<:WindowedArray, P<:Plan}

        nfft = length(plan)
        nout = div(nfft, 2) + 1

        win = winfun(winput.binsize)::Vector{Float64}

        norm2 = sum(abs2, win)
        r = fs * norm2
        @compat r_temp = Vector{Float64}(undef, 1)

        @compat fftbuf = Vector{Complex{Float64}}(undef, nout)
        winbuf = zeros(Float64, nfft)

        frequencies = collect(FFTW.rfftfreq(nfft, fs))

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
    winput::W, plan::P, fs::Real, winfun::Function = blackman
) where {W<:WindowedArray{<:Any, <:Any, 1, <:Any}, P<:Plan}
    Stft{Vector{Complex{Float64}}, W, P}(
        winput, plan, convert(Float64, fs), winfun
    )
end

function Stft(
    winput::W, fs::Real = 1, args...
) where W<:WindowedArray{<:Any, <:Any, 1, <:Any}
    nfft = nextfastfft(winput.binsize)
    @compat plan = plan_rfft(Vector{Float64}(undef, nfft))
    Stft(winput, plan, fs, args...)
end

function Stft(
    input::AbstractVector,
    binsize::Integer,
    fs::Real,
    overlap::Integer = 0,
    args...
)
    winput = WindowedArray(input, binsize, overlap)
    Stft(winput, fs, args...)
end

elsize(a::Stft) = (a.nout,)

frequencies(a::Stft) = a.frequencies

basedata(a::Stft) = basedata(a.winput)

function getindex(a::Stft, i::Integer)
    out = make_out(a)
    getindex!(out, a, i)
    out
end

function getindex!(out, a::Stft, i::Integer)
    @boundscheck checkbounds(a, i)
    window_index!(a, i)
    @static if VERSION >= v"0.7.0-DEV.2575"
        mul!(out, a.plan, a.winbuf)
    else
        A_mul_B!(out, a.plan, a.winbuf)
    end
    nothing
end

@compat make_out(a::Stft) = Vector{Complex{Float64}}(undef, a.nout)

function window_index!(a::Stft, i::Integer)
    v = a.winput[i]
    nv = length(v)

    if nv > 1
        if nv == length(a.win)
            win = a.win
            a.r_temp[1] = a.r
        else nv > 1
            win = a.winfun(nv)::Vector{Float64}
            a.r_temp[1] = a.fs * sum(abs2, win)
        end
        @simd for i in 1:nv
            @inbounds a.winbuf[i] = win[i] * v[i]
        end
    elseif nv == 1
        @inbounds a.winbuf[1] = 1
        a.r_temp[1] = a.fs
    else
        a.r_temp[1] = a.fs
    end
    @simd for i in (nv + 1):a.winput.binsize
        @inbounds a.winbuf[i] = 0
    end
end

function concrete_type(::Type{<:Stft}, ::Type{W}) where {W<:WindowedArray}
    return Stft{
        Vector{Float64},
        W,
        ConcretePlan
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
    stft::S
    fftbuf::Vector{Complex{Float64}}
    function StftPsd{E, S}(stft::S) where {E, S<:Stft}
        fftbuf = make_out(stft)
        new(stft, fftbuf)
    end
end
StftPsd(stft::S) where S<:Stft = StftPsd{Vector{Float64}, S}(stft)
StftPsd(args...) = StftPsd(Stft(args...))

@compat make_out(a::StftPsd) = Vector{Float64}(undef, a.stft.nout)

function getindex(s::StftPsd, i::Integer)
    dest = make_out(s)
    getindex_unsafe!(dest, s, i)
    dest
end

function getindex!(out, s::StftPsd, i::Integer)
    if ! copy_length_check(length(out), s.stft.nout)
        throw(ArgumentError("Dest must be at least as long as s_fft"))
    end
    getindex_unsafe!(out, s, i)
end

function copy!(
    out::AbstractArray,
    s::StftPsd,
    out_off::Integer = 1,
    s_off::Integer = 1,
    n_ndx::Integer = n_ndx(s_off, length(s))
)
    out_off > 0 && s_off > 0 || throw(ArgumentError("Invalid indices"))
    n_total = n_ndx * s.stft.nout
    rel_offset = (s_off - 1) * s.stft.nout + 1
    rel_len = length(s) * s.stft.nout
    if ! copy_length_check(length(out), rel_len, out_off, rel_offset, n_total)
        throw(ArgumentError("Destination is not large enough for copy!"))
    end
    cur_offset = out_off
    for i = 0:(n_ndx - 1)
        @inbounds getindex_unsafe!(out, s, s_off + i, cur_offset)
        cur_offset = cur_offset + s.stft.nout
    end
    out
end

function getindex_unsafe!(out, s::StftPsd, i::Integer, d_off::Integer = 1)
    getindex!(s.fftbuf, s.stft, i) # Does a bounds check
    fft2pow_unsafe!(
        out, s.fftbuf, s.stft.nfft, s.stft.r_temp[1], d_off, s.stft.nout
    )
end

frequencies(s::StftPsd) = frequencies(s.stft)
basedata(s::StftPsd) = basedata(s.stft)
elsize(a::StftPsd) = elsize(a.stft)

function fft2pow!(
    dest::AbstractVector,
    s_fft::AbstractVector,
    nfft::Integer,
    r::Real,
    d_offset::Integer = 1
)
    n = length(s_fft)
    if ! copy_length_check(length(dest), n)
        throw(ArgumentError("dest must be at least as long as s_fft"))
    end
    fft2pow_unsafe!(dest, s_fft, nfft, r, d_offset, n)
end

function fft2pow_unsafe!(
    dest::AbstractArray{T, <:Any},
    s_fft::AbstractVector{Complex{T}},
    nfft::Integer,
    r::Real,
    d_offset::Integer = 1,
    n::Integer = length(s_fft)
) where T
    m1 = convert(T, 1/r)
    m2 = 2 * m1
    @inbounds dest[d_offset] = m1 * abs2(s_fft[1])
    @simd for i = 1:(n - 2)
        @inbounds dest[d_offset + i] = m2 * abs2(s_fft[i])
    end
    @inbounds dest[d_offset + n - 1] = ifelse(iseven(nfft), m1, m2) * abs2(s_fft[n])
    nothing
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
