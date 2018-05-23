__precompile__()
module GLTimeseries

import Base: linearindexing, length, size, getindex, setindex!, extrema
import GLUtilities: bin_bounds, duration

using
    GLUtilities,
    DSP

export
    # Constants
    GL_CACHEPREFIX,
    # Types
    CachingDynamicTs,
    Downsampler,
    DynamicDownsampler,
    DynamicSpectrogram,
    DynamicTs,
    MappedDynamicDownsampler,
    MaxMin,

    # Functions
    cacheinfo,
    cache_reg,
    duration,
    downsamp_req,
    extent,
    make_shifter,
    open_cache_files,
    open_cache_file,
    parse_cache_filenames,
    shift_extrema,
    shift_extrema!,
    write_cache_files,
    write_cache_file

include("util.jl")
include("maxmin.jl")
include("dynamic_downsampler.jl")
include("dynamic_ts.jl")
include("mapped_ts.jl")
include("caching_ts.jl")
include("caching.jl")
include("spectrogram.jl")

end # module
