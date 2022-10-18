function u = get_exact(BLOCK, SUBDOMAIN, domain, OPTIONS)
    REAL = OPTIONS.REAL;
    x = BLOCK.grids.x;
    y = BLOCK.grids.y;
    e = SUBDOMAIN.edges;
    
    xt = @(x) x*(e(2)-e(1))/2 + (e(2)+e(1))/2;
    yt = @(y) y*(e(4)-e(3))/2 + (e(4)+e(3))/2;
    
    k = SUBDOMAIN.params.k(1);
    FUNC = OPTIONS.FUNC;
    
    u = zeros(length(x),length(y));
    
    for i = 1:length(x)
        for j = 1:length(y)
            if(domain(i,j))
                u(i,j) = u_exact(xt(x(i)), yt(y(j)), 0, 0, k, FUNC, REAL);
            end
        end
    end
    
end