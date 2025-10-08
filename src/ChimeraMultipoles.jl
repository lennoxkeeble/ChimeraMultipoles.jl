module ChimeraMultipoles
include("main.jl")

##### PRECOMPILATION TAKES 6000 SECS #####
using PrecompileTools: @setup_workload, @compile_workload    # this is a small dependency
using StaticArrays

# precompilation
@setup_workload begin
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
    d9x = @MArray [DX1(8), DX2(8), DX3(8)];

    q = 0.0042839;
    m = 1.0 + q;
    ν = q/(1+q)^2;

    i = 1; j = 2; k = 3; l = 1;
    OnePN, TwoPN, TwoPointFivePN = 1.0, 1.0, 1.0;
    

    @compile_workload begin
        SijDerivs.Sij2(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
        SijDerivs.Sij5(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
        SijDerivs.Sij6(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
        dxSijDerivs.dxkSij5(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)

        MijkDerivs.Mijk3(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
        MijkDerivs.Mijk7(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
        MijkDerivs.Mijk8(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
        dxMijkDerivs.dxlMijk7(i, j, k, l, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)

        MijDerivs.Mij2(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
        MijDerivs.Mij5(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
        MijDerivs.Mij6(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
        MijDerivs.Mij7(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
        MijDerivs.Mij8(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
        dxMijDerivs.dxkMij5(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
        dxMijDerivs.dxkMij6(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
        dxMijDerivs.dxkMij7(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)

    end
end
end