using BlindingIndex
using Test

@testset "BlindingIndex.jl" begin

    x = [48 22; 4 30; 330 319]

    @test addmargins(x) == [48 22 70; 4 30 34; 330 319 649; 382 371 753]
    @test nrow(x) == 3
    @test ncol(x) == 2
    @test bijames(x) == (bi_est = 0.896, bi_se = 0.011, bi_lcl = 0.875, bi_ucl = 0.918)
    @test bibang(x) == (bi_est = [0.115, 0.022], bi_se = [0.018, 0.019], bi_lcl = [0.08, -0.016], bi_ucl = [0.15, 0.06])
    
end
