module DynamicTimeseries

import Base: copy!, length, size, getindex, setindex!, extrema

import SignalIndices: bin_bounds, time_interval, duration

using DSP: DSP, blackman, nextfastfft, spectrogram, stft
using Distributed: Distributed, pmap
using FFTW: FFTW, AbstractFFTs, plan_rfft
using SignalIndices: SignalIndices, bin_center, clip_ndx, copy_length_check, div_type,
    n_ndx, ndx_to_t, t_to_last_ndx, t_to_ndx, view_trailing_slice, weighted_mean,
    weighted_mean_dim
using SortedIntervals: extrema_red
using LinearAlgebra: LinearAlgebra, mul!
using Mmap: Mmap
using EventIntervals: EventIntervals, MarkedPoint, NakedPoints, Point, Points,
    VariablePoints, point_values, points, pop_marks, pp_downsamp, pt_extent_merge,
    pt_merge
using Statistics: Statistics, mean
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


include("mmap_util.jl")
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
