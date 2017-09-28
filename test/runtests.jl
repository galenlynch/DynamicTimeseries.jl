using GLTimeseries
using Base.Test

@testset "GLTimeseries"  begin
    A = ones(1000)
    A[1:2:end] = 2

    @testset "timeseries" begin
        mm = MaxMin(A, 10)
        @test length(mm) == 100
        @test size(mm) == (100,)
        @test GLTimeseries.binsize(mm) == 10
        @test mm[1] == (1, 2)
        @test mm[end] == (1, 2)

        dts = DynamicTs(A, 10, 0)
        (xs, mm) = downsamp_req(dts, 0, 1, 10)
        (xs, mm) = downsamp_req(dts, 1000, 1001, 10) # Past signal
    end

    @testset "util" begin
        @test extent(A) == 1.0
        @test extent([A, A]) == [1.0, 1.0]
    end
end
