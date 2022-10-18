function GBf = get_GBf(BLOCK, SUBDOMAIN, REAL, ROD)
    
    if(ROD)
        GBf = get_GBf_rod(BLOCK, SUBDOMAIN, REAL);
    else
        GBf = get_GBf_square(BLOCK, SUBDOMAIN, REAL);
    end


end

function GBf = get_GBf_square(BLOCK, SUBDOMAIN, REAL)

    f = SUBDOMAIN.source.f;
    N1 = BLOCK.N(1); N2 = BLOCK.N(2);
    
    F = zeros(N1, N2);
    x = BLOCK.grids.x; y = BLOCK.grids.y; h = x(2) - x(1); 
    
    % Rather than loop through every element in both directions, we can
    % just loop through the x's and vectorize the evaluation of the y
    % values.
    for i = 1:N1
        for j = 1:N2
            F(i,j) = f(x(i),y(j),SUBDOMAIN.params.k(1));
        end
    end
    
    % Apply Bh and truncate to Mp
    Bf = apply_B(F,true);
    Bf(~BLOCK.grids.Mp) = 0;
    
    % Apply Gh and truncate to Np
    GBf = solver(Bf, SUBDOMAIN.params.k(1), BLOCK.N, h, true);
    GBf(~BLOCK.grids.Np) = 0;

end

function GBf = get_GBf_rod(BLOCK, SUBDOMAIN, REAL)

    f = SUBDOMAIN.source.f;
    N1 = BLOCK.N(1); N2 = BLOCK.N(2);
    F = zeros(N1, N2);
    x = BLOCK.grids.x; y = BLOCK.grids.y; h = x(2) - x(1);
    for i = 1:N1
        for j = 1:N2
            F(i,j) = f(x(i),y(j),SUBDOMAIN.params.k(2));
        end
    end
    Bf = apply_B(F,true);
    Bf(~BLOCK.grids.rodMp) = 0;
    GBf = solver(Bf, SUBDOMAIN.params.k(2), BLOCK.N, h, true);
    GBf(~BLOCK.grids.rodNp) = 0;
end


