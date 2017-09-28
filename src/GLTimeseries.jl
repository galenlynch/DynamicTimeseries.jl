module GLTimeseries

import Base: linearindexing, length, size, getindex, setindex!
import GLUtilities: bin_bounds

using GLUtilities

export
    # Types
    MaxMin,
    DynamicTs,
    Downsampler,

    # Functions
    extent,
    downsamp_req

include("util.jl")
include("timeseries.jl")

end # module
