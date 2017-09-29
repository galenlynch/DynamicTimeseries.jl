using GLTimeseries
using Base.Test

@testset "GLTimeseries"  begin
    const npt = 1000
    const A = ones(npt)
    A[1:2:end] = 2

    const fs = 10

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

        const maxpt = 10
        cdt = CachingDynamicTs(A, fs, 0, maxpt)
        (xs, mm) = downsamp_req(cdt, 0, 1, 100)
        @test length(mm) == 11
        (xs, mm) = downsamp_req(cdt, 0, 100, 100)
        @test length(mm) == 100
        (xs, mm) = downsamp_req(cdt, 0, 10, 50)
        @test length(mm) == 34
        (xs, mm) = downsamp_req(cdt, 0, 100, 10)
        @test length(mm) == 10
        cdt = CachingDynamicTs(A, fs, 0, maxpt, false)
        clean(cdt)
        cdt = 0
        cdt = CachingDynamicTs(A, fs, 0, maxpt, false)
        (cachepaths, cachelengths) = scavenge_cache(cdt)
        cdt2 = CachingDynamicTs(A, fs, 0, cachepaths, cachelengths)
        (xs, mm) = downsamp_req(cdt, 0, 1, 100)
        @test length(mm) == 11
    end

    @testset "spectrogram" begin
        ds = DynamicSpectrogram(A, 0, fs)
        (s, f, t) = downsamp_req(ds, 0, 100, 2)
    end

    @testset "util" begin
        @test extent(A) == 1.0
        @test extent([A, A]) == [1.0, 1.0]
    end
end
