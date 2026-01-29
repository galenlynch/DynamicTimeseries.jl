module GLTimeseries

import Base: copy!, length, size, getindex, setindex!, extrema, count

import GLUtilities: bin_bounds, time_interval, duration

using Compat, GLUtilities, DSP, PointProcesses, Missings

using FFTW, Statistics, Distributed, Mmap, LinearAlgebra
const Plan = AbstractFFTs.Plan
const ConcretePlan = FFTW.rFFTWPlan{Float64,-1,false,1}

export
    # Constants
    GL_CACHEPREFIX,
    # Types
    WindowedArray,
    DynamicWindower,
    CacheAccessor,
    CachingDynamicTs,
    CachingStftPsd,
    Downsampler,
    DynCachingStftPsd,
    AbstractDynamicDownsampler,
    DynamicDownsampler,
    AbstractDynamicSpectrogram,
    DynamicSpectrogram,
    MappedDynamicDownsampler,
    MaxMin,
    Averager,
    Stft,
    StftPsd,
    PointMerger,
    DynamicPointDownsampler,
    DynamicPointBoxer,

    # Functions
    basedata,
    cacheinfo,
    cache_reg,
    cachefiles_batch_mmap,
    duration,
    time_interval,
    downsamp_req,
    downsamp_batch_mmap,
    extent,
    frequencies,
    fs,
    getindex!,
    make_shifter,
    make_out,
    open_cache_files,
    open_cache_file,
    parse_cache_filenames,
    point_averager,
    shift_extrema,
    shift_extrema!,
    write_cache_files,
    write_cache_file


include("util.jl")
include("windower.jl")
include("abstract_downsamplers.jl")
include("dynamic_windower.jl")
include("dynamic_downsampler.jl")
include("maxmin.jl")
include("averager.jl")
include("stft.jl")
include("spectrogram.jl")
include("mapped_ts.jl")
include("cache_accessor.jl")
include("caching.jl")
include("caching_ts.jl")
include("caching_averager.jl")
include("caching_stft.jl")
include("dynamic_stft.jl")
include("batch.jl")
include("point_downsampler.jl")

end # module
