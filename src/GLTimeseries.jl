__precompile__()
module GLTimeseries

import Base: linearindexing, length, size, getindex, setindex!, extrema
import GLUtilities: bin_bounds, time_interval, duration, make_slice_idx

using
    GLUtilities,
    DSP

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
    AbstractDynamicDownsampler,
    DynamicDownsampler,
    DynamicSpectrogram,
    MappedDynamicDownsampler,
    MaxMin,
    Averager,
    Stft,
    StftPsd,

    # Functions
    cacheinfo,
    cache_reg,
    duration,
    time_interval,
    downsamp_req,
    extent,
    frequencies,
    getindex!,
    make_shifter,
    make_out,
    open_cache_files,
    open_cache_file,
    parse_cache_filenames,
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

end # module
