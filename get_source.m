function u = get_source(BLOCK, SUBDOMAIN, domain, OPTIONS)
    REAL = OPTIONS.REAL;
    x = BLOCK.grids.x;
    y = BLOCK.grids.y;
    k = SUBDOMAIN.params.k;
    FUNC = OPTIONS.FUNC;
    
    a = SUBDOMAIN.edges(1) - BLOCK.edges(1);
    b = SUBDOMAIN.edges(3) - BLOCK.edges(3);
    
    u = zeros(length(x),length(y));
    
    for i = 1:length(x)
        for j = 1:length(y)
            if(domain(i,j))
                if(SUBDOMAIN.type==2 && BLOCK.grids.rodMp(i,j))
                    k = SUBDOMAIN.params.krod;
                else
                    k = SUBDOMAIN.params.k;
                end
                u(i,j) = source(x(i), y(j), 0, 0, k, FUNC, REAL, false, a, b);
            end
        end
    end
    
end