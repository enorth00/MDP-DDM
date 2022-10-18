% set_OPTIONS.m

function [OPTIONS, BLOCK] = SET_GlobalOptions(OPTIONS)
    OPTIONS.BASIS = ["chebyshev", "chebyshev", "chebyshev", "chebyshev", "fourier"];
%     OPTIONS.Mmax = [101, 51];
%     OPTIONS.BF = set_basis_functions(OPTIONS.Mmax);
    OPTIONS.REAL = false;
    
    BLOCK = allocate_blocks(OPTIONS.types);
end