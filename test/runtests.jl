using GLTimeseries
using Base.Test

@testset "GLTimeseries"  begin

    const npt = 1000
    const A = ones(Int, npt)
    const B = rand(npt)
    const C = rand(Float32, npt)
    A[1:2:end] = 2

    const fs = 10

    @testset "timeseries" begin
        @testset "WindowedArray" begin
            wa = WindowedArray(A, 10)
            @test all(wa[1] == A[1:10])
            @test all(wa[2] == A[11:20])
            wa[end]
            wa_overlap = WindowedArray(A, 6, 1, 2)
            @test all(wa_overlap[1] == A[1:6])
            @test all(wa_overlap[2] == A[5:10])
            wab = WindowedArray(B, 10)
            @test all(wab[1] == B[1:10])
        end

        @testset "maxmin" begin
            mm = MaxMin(A, 10)
            @test length(mm) == 100
            @test size(mm) == (100,)
            @test GLTimeseries.binsize(mm) == 10
            @test mm[1] == (1, 2)
            @test mm[end] == (1, 2)
            mm[1:-1]
            mm2 = MaxMin(B, 10)
            mm2[1:-1]
            @test mm2[1] == extrema(B[1:10])
            arr2d = Array{Float64, 2}(2, length(mm2))
            for (i, (emin, emax)) in enumerate(mm2)
                arr2d[1, i] = emin
                arr2d[2, i] = emax
            end
            mm2d = MaxMin(arr2d, 10)
            mm2d[1:-1]
            real_e = (minimum(view(arr2d, 1, 1:10)), maximum(view(arr2d, 2, 1:10)))
            @test mm2d[1] == real_e
            mm_tup = MaxMin(mm2, 10)
            @test mm_tup[1] == real_e
        end

        @testset "averager" begin
            avg = Averager(B, 100)
            @test length(avg) == 10
            @test size(avg) == (10,)
            @test avg[1] == mean(B[1:100])
            @test avg[2] == mean(B[101:200])
            @test avg[1:2] == [mean(B[1:100]), mean(B[101:200])]
            avg_overlap = Averager(B, 7, 1, 3)
            @test avg_overlap[1] == mean(B[1:7])
            @test avg_overlap[2] == mean(B[5:11])
            D = rand(8, 4)
            println("D is ", D)
            avg_2d = Averager(D, 5, 1, 2)
            @test avg_2d[1] == mean(D[1:5, :], 1)
            @test avg_2d[2] == mean(D[4:8, :], 1)
            avg_2d_2 = Averager(D, 3, 2, 1)
            @test avg_2d_2[1] == mean(D[:, 1:3], 2)
            avg_int = Averager(A, 10)
            @test avg_int[1] == mean(A[1:10])
        end

        @testset "DynamicWindower" begin
            dw = DynamicWindower(B, fs)
            (xs, ys, wd) = downsamp_req(dw, 0, 10, 2)
            (xs, ys, wd) = downsamp_req(dw, 0, 10, 0)
        end
        @testset "DynamicTs" begin
            dts = DynamicTs(A, 10, 0) # 10 Hz 0 offset, should be 100 s of signal
            (xs, mm, was_downsamped) = downsamp_req(dts, 1, 2, 10) # Get 10 values from 0 - 1 seconds
            @test length(mm) == 7
            (xs, mm, was_downsamped) = downsamp_req(dts, 1000, 1001, 10) # Past signal
            @test length(mm) == 0
            (xs, mm, wd) = downsamp_req(dts, 0, 0, 1)
            (xs, mm, wd) = downsamp_req(dts, 0, 1, 0)
        end

        @testset "CachingDynamicTs" begin
            const maxpt = 10
            @testset "Caching" begin
                cachetype = CachingDynamicTs{Int, typeof(A)}
                (path, npair) = write_cache_file(
                    cachetype, A
                )
                (files, lengths) = write_cache_files(
                    cachetype, A, 10
                )
                println(files)
                println(lengths)
                cachedir = tempdir()
                (files, lengths) = write_cache_files(
                    cachetype, A, 10, true;
                    cachedir=cachedir, fid=1)
                println(files)
                println(lengths)
                cachearr = open_cache_files(cachetype, cachedir, 1)
            end

            cdt = CachingDynamicTs(A, fs, 0, maxpt)
            (xs, mm, wd) = downsamp_req(cdt, 0, 0, 1)
            (xs, mm, wd) = downsamp_req(cdt, 0, 1, 0)
            (xs, mm, was_downsamped) = downsamp_req(cdt, 0, 1, 100)
            @test length(mm) == 11
            println(time_interval(cdt))
            (xs, mm, wd) = downsamp_req(cdt, 0, 100, 10, false)
            (xs, mm, was_downsamped) = downsamp_req(cdt, 0, 100, 100)
            @test length(mm) == 100
            (xs, mm, was_downsamped) = downsamp_req(cdt, 0, 10, 50)
            @test length(mm) == 51
            (xs, mm, was_downsamped) = downsamp_req(cdt, 0, 100, 10)
            @test length(mm) == 10
            cdt = CachingDynamicTs(A, fs, 0, maxpt, false)
            const npt_C = 10000
            const C = rand(10000)
            const fs_C = 100
            const dts_C = CachingDynamicTs(C, fs_C)
            const true_extrema = extrema(C)
            (xs, ys, _) = downsamp_req(dts_C, 0.0, 99.99, 74)
            @test GLTimeseries.extrema_red(ys) == true_extrema
            @testset "MappedDynamicDownsampler" begin
                c = 5
                mds = MappedDynamicDownsampler(cdt, make_shifter(c))
                (xs, ys) = downsamp_req(mds, 0, 100, 10)
                (xs, mm, wd) = downsamp_req(mds, 0, 0, 1)
                (xs, mm, wd) = downsamp_req(mds, 0, 1, 0)
                println(mm)
                println(ys)
            end
        end
    end

    @testset "spectrogram" begin
        const ds = DynamicSpectrogram(A, fs, 0)
        (t, (f, s, t_w, f_w)) = downsamp_req(ds, 0, 100, 2)
        (t, (f, s, t_w, f_w)) = downsamp_req(ds, 0, 0, 2)
        (t, (f, s, t_w, f_w)) = downsamp_req(ds, 0, 100, 0)
        const ds32 = DynamicSpectrogram(C, fs, 0)
        (t, (f, s, t_w, f_w)) = downsamp_req(ds32, 0, 100, 0)
    end

    @testset "util" begin
        @test extent(A) == 1.0
        @test extent([A, A]) == [1.0, 1.0]
    end
end
