using GLTimeseries
using Base.Test

@testset "GLTimeseries"  begin
    @testset "timeseries" begin
        A = rand(1000, 1)
        dts = DynamicTs(A, 10, 0)
        (xs, mm) = downsamp_req(dts, 1, 2, 20)
    end
end
