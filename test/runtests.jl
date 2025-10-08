using ChimeraMultipoles
using Test
using StaticArrays

## TEST RADIATION REACTION POTENTIALS ##
DX1(n) = sin(n+1) * cos(n+1)
DX2(n) = sin(5n+2) * cos(5n+2)
DX3(n) = sin(10n+3) * cos(10n+3)

x = @MArray [DX1(0), DX2(0), DX3(0)];
dx = @MArray [DX1(1), DX2(1), DX3(1)];
d2x = @MArray [DX1(2), DX2(2), DX3(2)];
d3x = @MArray [DX1(3), DX2(3), DX3(3)];
d4x = @MArray [DX1(4), DX2(4), DX3(4)];
d5x = @MArray [DX1(5), DX2(5), DX3(5)];
d6x = @MArray [DX1(6), DX2(6), DX3(6)];
d7x = @MArray [DX1(7), DX2(7), DX3(7)];
d8x = @MArray [DX1(8), DX2(8), DX3(8)];
d9x = @MArray [DX1(9), DX2(9), DX3(9)];

q = 0.0042839;
m = 1.0 + q;
ν = q/(1+q)^2;

i = 1; j = 2; k = 3; l = 1;
OnePN, TwoPN, TwoPointFivePN = 1.0, 1.0, 1.0;

@time ChimeraMultipoles.SijDerivs.Sij2(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
@time ChimeraMultipoles.SijDerivs.Sij5(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
@time ChimeraMultipoles.SijDerivs.Sij6(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
@time ChimeraMultipoles.dxSijDerivs.dxkSij5(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)

@time ChimeraMultipoles.MijkDerivs.Mijk3(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
@time ChimeraMultipoles.MijkDerivs.Mijk7(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
@time ChimeraMultipoles.MijkDerivs.Mijk8(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
@time ChimeraMultipoles.dxMijkDerivs.dxlMijk7(i, j, k, l, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)

@time ChimeraMultipoles.MijDerivs.Mij2(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
@time ChimeraMultipoles.MijDerivs.Mij5(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
@time ChimeraMultipoles.MijDerivs.Mij6(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
@time ChimeraMultipoles.MijDerivs.Mij7(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
@time ChimeraMultipoles.MijDerivs.Mij8(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
@time ChimeraMultipoles.dxMijDerivs.dxkMij5(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
@time ChimeraMultipoles.dxMijDerivs.dxkMij6(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
@time ChimeraMultipoles.dxMijDerivs.dxkMij7(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)

MijDerivs_MMA = [-1.2731132054647223, -146.22391259001472, -1047.3312753529622, -7269.584662808288, -24206.498471952626];
dxMijDerivs_MMA = [-838.2900612190671, -1286.864884380358, 58728.78156599056];
MijkDerivs_MMA = [-1.5729103024914632, 3417.04201941972, 36168.57745813809];
dxMijkDeriv = -21297.57952429885;
SijDerivs_MMA = [0.4327051584775146, -19.7993591714404, -71.78753092185924];
dxSijDeriv = -82.40911402192786;

MijDerivs_julia = [ChimeraMultipoles.MijDerivs.Mij2(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν), 
    ChimeraMultipoles.MijDerivs.Mij5(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν), 
    ChimeraMultipoles.MijDerivs.Mij6(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν), 
    ChimeraMultipoles.MijDerivs.Mij7(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν), 
    ChimeraMultipoles.MijDerivs.Mij8(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)];

dxMijDerivs_julia = [ChimeraMultipoles.dxMijDerivs.dxkMij5(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν),
    ChimeraMultipoles.dxMijDerivs.dxkMij6(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν),
    ChimeraMultipoles.dxMijDerivs.dxkMij7(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)];


SijDerivs_julia = [ChimeraMultipoles.SijDerivs.Sij2(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν),
    ChimeraMultipoles.SijDerivs.Sij5(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν),
    ChimeraMultipoles.SijDerivs.Sij6(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)];

dxSijDeriv_julia = ChimeraMultipoles.dxSijDerivs.dxkSij5(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν);

MijkDerivs_julia = [ChimeraMultipoles.MijkDerivs.Mijk3(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν),
    ChimeraMultipoles.MijkDerivs.Mijk7(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν),
    ChimeraMultipoles.MijkDerivs.Mijk8(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)];

dxMijkDeriv_julia = ChimeraMultipoles.dxMijkDerivs.dxlMijk7(i, j, k, l, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν);


@testset "ChimeraMultipoles.jl" begin
    @test isapprox(MijDerivs_MMA, MijDerivs_julia, rtol = 1e-12)
    @test isapprox(dxMijDerivs_MMA, dxMijDerivs_julia, rtol = 1e-12)
    @test isapprox(SijDerivs_MMA, SijDerivs_julia, rtol = 1e-12)
    @test isapprox(dxSijDeriv, dxSijDeriv_julia, rtol = 1e-12)
    @test isapprox(MijkDerivs_MMA, MijkDerivs_julia, rtol = 1e-12)
    @test isapprox(dxMijkDeriv, dxMijkDeriv_julia, rtol = 1e-12)
end