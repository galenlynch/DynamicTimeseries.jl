abstract type PointDownsampler{E} <: AbstractDynamicDownsampler{E} end

struct MarkedPointMerger{A<:MarkedPoint,R,P<:MarkedPoint{<:Any,R}}
    mergefunc::Function
    outmap::Function
end

function point_extent_averager(::Type{P}) where {E<:Real,M<:Real,P<:MarkedPoint{E,M}}
    F = div_type(E)
    L = div_type(M)
    R = Tuple{Tuple{E,E},L}
    MarkedPointMerger{P,R,MarkedPoint{F,R}}(pt_extent_merge, pop_marks)
end

function point_averager(::Type{P}) where {E<:Real,M<:Real,P<:MarkedPoint{E,M}}
    F = div_type(E)
    R = div_type(M)
    MarkedPointMerger{P,R,MarkedPoint{F,R}}(pt_merge, pop_marks)
end
point_averager(::P) where {P<:Point} = point_averager(P)

# For some reason the compilter can't reason about this output, even though
# the function call itself is fine
function merge_points(
    pm::MarkedPointMerger{A,<:Any,P},
    pts::Points{<:Any,1,<:Any,A},
    xb,
    xe,
    resolution,
) where {A,P}
    pm.outmap(pp_downsamp(pts, xb, xe, resolution, pm.mergefunc, P))
end

# Version to help the compiler do type inference
function merge_points(
    pm::MarkedPointMerger{A,<:Any,P},
    pts::VariablePoints{E,I,<:Any,M,<:Any},
    xb,
    xe,
    resolution,
) where {E,I,M,A<:MarkedPoint{E,M},F,N,P<:MarkedPoint{F,N}}
    pm.outmap(
        pp_downsamp(pts, xb, xe, resolution, pm.mergefunc, P),
    )::VariablePoints{F,I,NakedPoints{F,I,Vector{F}},N,Vector{N}}
end


struct DynamicPointDownsampler{
    PointDimension,
    PointType<:MarkedPoint{PointDimension,<:Any},
    MergedMarkType,
    MergerType<:MarkedPointMerger{PointType,MergedMarkType,<:Any},
    P<:Points{PointDimension,1,<:Any,PointType},
} <: PointDownsampler{MergedMarkType}
    points::P
    merger::MergerType
end

@static if VERSION < v"0.7.0-DEV.2575"
    function DynamicPointDownsampler(
        points::P,
        merger::M,
    ) where {
        F,
        A<:MarkedPoint{F,<:Any},
        E,
        M<:MarkedPointMerger{A,E,<:Any},
        P<:Points{F,1,<:Any,A},
    }
        DynamicPointDownsampler{F,A,E,M,P}(points, merger)
    end
end

function DynamicPointDownsampler(points::Points{<:Number,1,<:Any,M}) where {M<:MarkedPoint}
    DynamicPointDownsampler(points, point_extent_averager(M))
end

basepoints(dp::DynamicPointDownsampler) = dp.points

fs(::DynamicPointDownsampler) = error("Not Implemented")
baselength(::DynamicPointDownsampler) = error("Not Implemented")

start_time(dp::DynamicPointDownsampler) = time_interval(dp.points)[1]
time_interval(dp::DynamicPointDownsampler) = time_interval(dp.points)
duration(dp::DynamicPointDownsampler) = duration(dp.points)

function downsamp_req(dp::DynamicPointDownsampler, xb, xe, resolution)
    merged = merge_points(dp.merger, dp.points, xb, xe, resolution)
    (point_values(merged)..., true)
end

struct DynamicPointBoxer{
    PtDomainT,
    PtRangeT,
    PtType<:Tuple{<:NTuple{2,PtDomainT},PtRangeT},
    D<:DynamicPointDownsampler{<:Any,<:Any,PtType,<:Any,<:Any},
} <: PointDownsampler{PtType}
    pointdownsampler::D
    min_width::PtDomainT
    y_center::PtRangeT
end

function DynamicPointBoxer(
    dp::D,
    min_width::Number,
    y_center::Number = 0,
) where {
    PtDomainT,
    PtRangeT,
    PtType<:Tuple{<:NTuple{2,PtDomainT},PtRangeT},
    D<:DynamicPointDownsampler{<:Any,<:Any,PtType},
}
    DynamicPointBoxer{PtDomainT,PtRangeT,PtType,D}(
        dp,
        convert(PtDomainT, min_width),
        convert(PtRangeT, y_center),
    )
end

function DynamicPointBoxer(
    points::Points{<:Number,1,<:Any,M},
    args...,
) where {M<:MarkedPoint}
    DynamicPointBoxer(DynamicPointDownsampler(points, point_extent_averager(M)), args...)
end

basepoints(dp::DynamicPointBoxer) = basepoints(dp.pointdownsampler)

function downsamp_req(
    dp::DynamicPointBoxer{PtDomainT,PtRangeT,<:Any,<:Any},
    xb,
    xe,
    res,
) where {PtDomainT,PtRangeT}
    # ptvals will start as center points for the box,
    # but I will make it the left point inplace
    ptvals, marks, _ = downsamp_req(dp.pointdownsampler, xb, xe, res)
    npt = length(ptvals)
    @compat heights = Vector{PtRangeT}(undef, npt)
    bottoms = similar(heights)
    @compat widths = Vector{PtDomainT}(undef, npt)
    min_half_width = dp.min_width / 2
    for (i, ((bb, be), height)) in enumerate(marks)
        bottoms[i] = dp.y_center - height / 2
        heights[i] = height
        raw_width = be - bb
        if raw_width >= dp.min_width
            ptvals[i] = bb - min_half_width # Set this to the 'left' part of the range
            widths[i] = raw_width + dp.min_width
        else
            ptvals[i] = ptvals[i] - min_half_width
            widths[i] = dp.min_width
        end
    end
    ptvals, (bottoms, widths, heights), true
end

fs(::DynamicPointBoxer) = error("Not Implemented")
baselength(::DynamicPointBoxer) = error("Not Implemented")

function extrema(dp::DynamicPointBoxer{<:Any,PtRangeT,<:Any,<:Any}) where {PtRangeT}
    _, heights = point_values(basepoints(dp))
    if isempty(heights)
        res = (dp.y_center, dp.y_center)
    else
        half_height = maximum(heights) / 2
        res = (dp.y_center - half_height, dp.y_center + half_height)
    end
    res
end

start_time(dp::DynamicPointBoxer) = start_time(dp.pointdownsampler)
time_interval(dp::DynamicPointBoxer) = time_interval(dp.pointdownsampler)
duration(dp::DynamicPointBoxer) = duration(dp.pointdownsampler)
