%-------------------------------------------------------------------------
%
%   Function: set_basis_functions
%
%   Inputs:     MMax - contains [McMax, MrMax] which are typically the
%                       maximum number of basis functions we could need
%                       for the sides (McMax) or rod (MrMax).
%               types - two arguments, [sides, rod]. Should be 'chebyshev'
%                       or 'fourier'.
%
%   Outputs:   bf - object containing the precomputed matrices of data
%
%   Description:    This function precomputes handles for the basis 
%                   functions to be used for the problem. This includes 
%                   the second and fourth derivatives.
%
%-------------------------------------------------------------------------


function bf = SET_basis_functions(MMax, types, OPTIONS)
    Mc = MMax(1);
    Mr = MMax(2);
    p = OPTIONS.FOURIER_AUX;
    % Pick the appropriate basis functions
    if(types == "chebyshev")
        b0 = chebpoly(0:(Mc-1), [-1,1]);
    elseif(types == "fourier")
        b0 = trigpoly((Mc-1)/2:-1:(-(Mc-1)/2), p);
    end
    
    % Regardless of the type (cheb vs. trig), diff can handle the rest.
    b2 = diff(b0,2);
    b4 = diff(b2,2);

    bf.b0 = b0;
    bf.b2 = b2;
    bf.b4 = b4;

    % These are for the rod. Always Fourier (for now at least)
    f0 = trigpoly((Mr-1)/2:-1:(-(Mr-1)/2), [-pi, pi]);
    f2 = diff(f0, 2);
    f4 = diff(f2, 2);
    
    r.d0.xi0 = f0;
    r.d0.xi1 = f0;
    r.d2.xi0 = f2;
    r.d2.xi1 = f2;
    r.d4.xi0 = f4;
    
    bf.r = r;
end
