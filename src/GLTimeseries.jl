module GLTimeseries

using GLUtilities

import Base: linearindexing, length, size, getindex, setindex!

export
    # Types
    MaxMin,
    DynamicTs,
    Downsampler,

    # Functions
    downsamp_req

include("timeseries.jl")

end # module
