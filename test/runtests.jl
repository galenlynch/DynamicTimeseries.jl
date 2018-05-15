using GLTimeseries
using Base.Test

@testset "GLTimeseries"  begin
    const npt = 1000
    const A = ones(Int, npt)
    const B = rand(npt)
    A[1:2:end] = 2

    const fs = 10

    @testset "timeseries" begin
        @testset "maxmin" begin
            mm = MaxMin(A, 10)
            @test length(mm) == 100
            @test size(mm) == (100,)
            @test GLTimeseries.binsize(mm) == 10
            @test mm[1] == (1, 2)
            @test mm[end] == (1, 2)
            mm2 = MaxMin(B, 10)
            @test mm2[1] == extrema(B[1:10])
            arr2d = Array{Float64, 2}(2, length(mm2))
            for (i, (emin, emax)) in enumerate(mm2)
                arr2d[1, i] = emin
                arr2d[2, i] = emax
            end
            mm2d = MaxMin(arr2d, 10)
            real_e = (minimum(view(arr2d, 1, 1:10)), maximum(view(arr2d, 2, 1:10)))
            @test mm2d[1] == real_e
            mm_tup = MaxMin(mm2, 10)
            @test mm_tup[1] == real_e
        end

        dts = DynamicTs(A, 10, 0)
        (xs, mm, was_downsamped) = downsamp_req(dts, 0, 1, 10)
        (xs, mm, was_downsamped) = downsamp_req(dts, 1000, 1001, 10) # Past signal

        const maxpt = 10
        cdt = CachingDynamicTs(A, fs, 0, maxpt)
        (xs, mm, was_downsamped) = downsamp_req(cdt, 0, 1, 100)
        @test length(mm) == 11
        (xs, mm, was_downsamped) = downsamp_req(cdt, 0, 100, 100)
        @test length(mm) == 100
        (xs, mm, was_downsamped) = downsamp_req(cdt, 0, 10, 50)
        @test length(mm) == 34
        (xs, mm, was_downsamped) = downsamp_req(cdt, 0, 100, 10)
        @test length(mm) == 10
        cdt = CachingDynamicTs(A, fs, 0, maxpt, false)
        (path, npair) = write_cache_file(A)
        (files, lengths) = write_cache_files(A, 10)
        println(files)
        println(lengths)
        cachedir = tempdir()
        (files, lengths) = write_cache_files(cachedir, 1, true, A, 10)
        println(files)
        println(lengths)
        cachearr = open_cache_files(Int, cachedir, 1)

        c = 5
        mds = MappedDynamicDownsampler(cdt, make_shifter(c))
        (xs, ys) = downsamp_req(mds, 0, 100, 10)
        println(mm)
        println(ys)
    end

    @testset "spectrogram" begin
        ds = DynamicSpectrogram(A, fs, 0)
        (t, (f, s)) = downsamp_req(ds, 0, 100, 2)
    end

    @testset "util" begin
        @test extent(A) == 1.0
        @test extent([A, A]) == [1.0, 1.0]
    end
end
