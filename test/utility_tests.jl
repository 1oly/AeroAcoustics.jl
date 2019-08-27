using Test
using AeroAcoustics

oct1 = [(11,22),
    (22,44),
    (44,88),
    (88,177),
    (177,355),
    (355,710),
    (710,1420),
    (1420,2840),
    (2840,5680),
    (5680,11360),
    (11360,22720)]


@testset "SPL:" begin
    @test isnan(SPL(-1.0))
    @test isnan(SPL(NaN))
    @test isnan(SPL(missing))
    @test SPL(1.0) ≈ 94 rtol=1e-3
    @test SPL.([1,4]) ≈ [94,100] rtol=1e-3
end

@testset "octavebands:" begin
    # Test 1/1 octavebands (using oct1 table above)
    fc = AeroAcoustics.octavebands_nomial(1)
    d = octavebandlimits(fc,1)
    dt = collect(zip(d...))
    dtInt = map(x->round.(Int,x), dt)
    for i in 1:length(oct1)
        @test all(isapprox.(oct1[i],dtInt[i];rtol=0.05))
    end
end

@testset "narrow2oct:" begin
    n = 100
    fc = octavebands(1;nomial=true)
    x = FreqArray(ones(n,length(fc)),fc)
    @test narrow2oct(x,1).arr == ones(n,length(fc))
    @test narrow2oct(x,1).fc == fc
end

@testset "point_to_region:" begin
    pt = (0.0,0.0)
    @test AeroAcoustics.point_to_region(pt,(1.0,2.0)) == [-0.5,0.5,-1.0,1.0]
    @test AeroAcoustics.point_to_region(pt,(1.0,)) == [-0.5,0.5,-0.5,0.5]
    @test AeroAcoustics.point_to_region(pt,1.0) == [-0.5,0.5,-0.5,0.5]
    @test AeroAcoustics.point_to_region(pt,1) == [-0.5,0.5,-0.5,0.5]
end
