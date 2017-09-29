module GLTimeseries

import Base: linearindexing, length, size, getindex, setindex!
import GLUtilities: bin_bounds

using GLUtilities

export
    # Types
    MaxMin,
    DynamicTs,
    DynamicDownsampler,
    CachingDynamicTs,
    Downsampler,

    # Functions
    extent,
    downsamp_req,
    clean,
    scavenge_cache

include("util.jl")
include("timeseries.jl")

end # module
