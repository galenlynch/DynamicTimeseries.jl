module GLTimeseries

import Base: linearindexing, length, size, getindex, setindex!
import GLUtilities: bin_bounds

using
    GLUtilities,
    DSP

export
    # Types
    MaxMin,
    DynamicTs,
    DynamicDownsampler,
    DynamicSpectrogram,
    CachingDynamicTs,
    Downsampler,

    # Functions
    extent,
    downsamp_req,
    clean,
    scavenge_cache

include("util.jl")
include("timeseries.jl")
include("spectrogram.jl")

end # module
