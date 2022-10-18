function [F, GBf, F_ext, F_rod, rod_ext] = source_extender(BLOCK, SUBDOMAIN, OPTIONS)
    REAL = OPTIONS.REAL;

    F_ext = extend_xi_inhomogeneous(BLOCK, SUBDOMAIN);
    F_ext_Matrix = zeros(BLOCK.N(1), BLOCK.N(2));
    F_ext_Matrix(BLOCK.grids.Ngamma) = F_ext;

    [A_f, ~] = N_plus_potential(F_ext_Matrix, BLOCK, SUBDOMAIN, REAL);
    
    GBf = get_GBf(BLOCK, SUBDOMAIN, REAL, false);
    
    F = F_ext - GBf(BLOCK.grids.Ngamma) - A_f(BLOCK.grids.Ngamma);
    
    if(SUBDOMAIN.params.M(5)>0)
        rod_ext = extend_xi_inhomogeneous_rod(BLOCK,SUBDOMAIN);
        rod_matrix = zeros(size(F_ext_Matrix));
        rod_matrix(BLOCK.grids.rodNgamma) = rod_ext;
        [A_rod, ~] = N_plus_potential_rod(rod_matrix, BLOCK, SUBDOMAIN, REAL);
        GBfrod = get_GBf(BLOCK, SUBDOMAIN, REAL, true);
        
        F_rod = rod_ext - GBfrod(BLOCK.grids.rodNgamma) - A_rod(BLOCK.grids.rodNgamma);
    else
        F_rod = [];
        rod_ext = [];
    end
    
    
end