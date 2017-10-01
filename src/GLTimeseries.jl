__precompile__()
module GLTimeseries

import Base: linearindexing, length, size, getindex, setindex!
import GLUtilities: bin_bounds, duration

using
    GLUtilities,
    DSP

export
    # Constants
    GL_CACHEPREFIX,
    # Types
    MaxMin,
    DynamicTs,
    DynamicDownsampler,
    DynamicSpectrogram,
    CachingDynamicTs,
    MappedDynamicDownsampler,
    Downsampler,

    # Functions
    extent,
    duration,
    downsamp_req,
    write_cache_files,
    write_cache_file,
    open_cache_files,
    open_cache_file,
    shift_extrema,
    shift_extrema!,
    make_shifter,
    cache_reg

include("util.jl")
include("timeseries.jl")
include("spectrogram.jl")
include("caching.jl")

end # module
