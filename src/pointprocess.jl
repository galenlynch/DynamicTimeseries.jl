# If not SimplePointProcess, expected to have a field "pointprocess" that is a
# SimplePointProcess or implement a pointprocess method that will yield a
# simplepointprocess
# Alternatively, implement duration, point_range, and count
abstract type PointProcess{T, E} end

rate(p::PointProcess) = count(p) / duration(p)
rate(p::PointProcess, b, e) = count(p, b, e) / (e - b)

struct SimplePointProcess{E<:Number, A<:AbstractVector{E}} <: PointProcess{E, Void}
    interval::NTuple{2, E}
    points::A
    function SimplePointProcess{E,A}(
        points::A, interval::NTuple{2, E}
    ) where {E, A<:AbstractVector{E}}
        issorted(points) || throw(ArgumentError("Times is not sorted"))
        validate_interval(interval) || throw(ArgumentError("Interval not valid"))
        if ! isempty(points)
            if points[1] < interval[1] || points[end] > interval[2]
                throw(ArgumentError("Points exceed interval"))
            end
        end
        new(interval, points)
    end
end

function SimplePointProcess(
    points::A, interval::NTuple{2, E}, check::Bool = true
) where {E, A<:AbstractVector{E}}
    if check && ! issorted(points)
        sort!(points)
    end
    SimplePointProcess{E, A}(points, interval)
end

function SimplePointProcess(
    points::AbstractVector{E}, interval::NTuple{2, <:Any}
) where {E}
    SimplePointProcess(points, (convert.(E, interval)...))
end

function SimplePointProcess(points::AbstractVector{E}, a::Number, b::Number) where E
    SimplePointProcess(points, (convert(E, a), convert(E, b)))
end

function SimplePointProcess(points::AbstractVector)
    isempty(points) && throw(ArgumnetError("points cannot be empty"))
    sort!(points)
    SimplePointProcess(points, (points[1], points[end]), false)
end

validate_interval(int::NTuple{2, Number}) = int[2] >= int[1]

count(p::SimplePointProcess) = length(p.points)
function count(p::SimplePointProcess, range_start, range_end)
    ib, ie = point_range(p, range_start, range_end)
    max(zero(ib), n_ndx(ib, ie))
end

function point_range(p::SimplePointProcess{E, <:Any}, b::E, e::E) where E
    b <= e || throw(ArgumentError("beginning and end not well ordered"))
    ib = searchsortedfirst(p.points, b)
    ie = searchsortedlast(p.points, e)
    ib, ie
end

duration(p::SimplePointProcess) = p.interval[2] - p.interval[1]

function points(p::SimplePointProcess, b, e)
    ib, ie = point_range(p, b, e)
    view(p.points, ib:ie)
end
points(p::SimplePointProcess) = p.points

pointprocess(p::SimplePointProcess) = p
pointprocess(p::PointProcess) = p.pointprocess

duration(p::PointProcess) = duration(pointprocess(p))
point_range(p::PointProcess, args...) = point_range(pointprocess(p), args...)
count(p::PointProcess, args...) = count(pointprocess(p), args...)

abstract type MarkedPointProcess{E, M} <: PointProcess{E, M} end

struct VariablePointProcess{
    E, P<:SimplePointProcess{E}, M, A<:AbstractVector{M}
} <: MarkedPointProcess{E, M}
    pointprocess::P
    marks::A
    function VariablePointProcess{E, P, M, A}(
        pointprocess::P, marks::A
    ) where {E, P<:SimplePointProcess{E}, M, A<:AbstractVector{M}}
        if count(pointprocess) != length(marks)
            throw(DimensionMismatch("point process must be the same size as marks"))
        end
        new(pointprocess, marks)
    end
end

function VariablePointProcess(
    pointprocess::P, marks::A
) where {E, P<:SimplePointProcess{E}, M, A<:AbstractVector{M}}
    VariablePointProcess{E, P, M, A}(pointprocess, marks)
end

function VariablePointProcess(
    points::AbstractVector, marks::AbstractVector, args...
)
    pp = SimplePointProcess(points, args...)
    VariablePointProcess(pp, marks)
end

function points(mp::MarkedPointProcess, ib, ie)
    pp = pointprocess(mp)
    ib, ie = point_range(pp, ib, ie)
    view(pp.points, ib:ie), view(mp.marks, ib:ie)
end
points(mp::MarkedPointProcess) = points(mp.pointprocess), mp.marks

struct SubPointProcess{E, M, P<:PointProcess{E, M}} <: PointProcess{E, M}
    pointprocess::P
    interval::NTuple{2, E}
    function SubPointProcess{E, M, P}(
        pointprocess::P, interval::NTuple{2, E}
    ) where {E, M, P<:PointProcess{E, M}}
        if ! validate_interval(interval)
            throw(ArgumentError("Invalid interval"))
        end
        if interval[1] < pointprocess.interval[1] ||
            interval[2] > pointprocess.interval[2]
            throw(ArgumentError("sub interval is not contained in parent interval"))
        end
        new(pointprocess, interval)
    end
end

function SubPointProcess(
    pointprocess::P, interval::NTuple{2, E}
) where {E, M, P<:PointProcess{E, M}}
    SubPointProcess{E, M, P}(pointprocess, interval)
end

function SubPointProcess(
    pointprocess::PointProcess{E, <:Any}, interval::NTuple{2, <:Any}
) where E
    SubPointProcess(pointprocess, (convert.(E, interval)...))
end

function SubPointProcess(pointprocess::PointProcess{E, <:Any}, b, e) where E
    SubPointProcess(pointprocess, (convert(E, b), convert(E, e)))
end

duration(spp::SubPointProcess) = spp.interval[2] - spp.interval[1]
count(spp::SubPointProcess) = count(spp.pointprocess, spp.interval...)
function count(spp::SubPointProcess, b, e)
    int_int = interval_intersect(spp.interval..., b, e)
    isempty(int_int) ? zero(b) : count(spp.pointprocess, int_int...)
end
    
function points(spp::SubPointProcess, b, e)
    int_int = interval_intersect(spp.interval..., b, e)
    if isempty(int_int)
        res = points(spp.pointprocess, 1, 0)
    else
        res = points(spp.pointprocess, int_int...)
    end
    res
end
