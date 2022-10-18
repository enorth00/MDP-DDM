
% Debugging purposes.
% 
% Expands the exact (known) solution and flux, passing along the
% coefficients instead of solving for the coefficients in the system.

function c = DEBUG_get_exact_coeffs(BLOCK, SUBDOMAIN, OPTIONS)

% debugging purposes, evaluates the exact coefficients for dirichlet and
% neumann data.

    FUNC = OPTIONS.FUNC;
    l = SUBDOMAIN.edges;
    k = SUBDOMAIN.params.k(1);
    M = SUBDOMAIN.params.M(1:4);
    Mr = SUBDOMAIN.params.M(5);
    p = pi/2 - 1;
    P = [l(1)-p, l(2)+p, l(3)-p, l(4)+p];
    
    a1 = BLOCK.edges(1); c1 = BLOCK.edges(3);
    xt = @(x) x + (l(1) - a1);
    yt = @(y) y + (l(3) - c1);
    
    rad = BLOCK.radius;
    REAL = OPTIONS.REAL;
    
    
    if(SUBDOMAIN.basis(1) == 'c')
        c01 = chebcoeffs(chebfun(@(x) u_exact(l(1),yt(x),0,0,k(1),FUNC,REAL,false), [-1,1]), M(1));
        c11 = chebcoeffs(chebfun(@(x) -u_exact(l(1),yt(x),1,0,k(1),FUNC,REAL,false), [-1,1]), M(1));
    elseif(SUBDOMAIN.basis(1) == 'f')
        c01 = trigcoeffs(chebfun(@(x) u_exact(l(1),yt(x),0,0,k(1),FUNC,REAL,false), [-pi/2, pi/2], 'trig'), M(1));
        c11 = trigcoeffs(chebfun(@(x) -u_exact(l(1),yt(x),1,0,k(1),FUNC,REAL,false), [-pi/2, pi/2], 'trig'), M(1));
    end

    if(SUBDOMAIN.basis(2) == 'c')
        c02 = chebcoeffs(chebfun(@(x) u_exact(l(2),yt(x),0,0,k(1),FUNC,REAL,false), [-1,1]), M(2));
        c12 = chebcoeffs(chebfun(@(x) u_exact(l(2),yt(x),1,0,k(1),FUNC,REAL,false), [-1,1]), M(2));
    elseif(SUBDOMAIN.basis(2) == 'f')
        c02 = trigcoeffs(chebfun(@(x) u_exact(l(2),yt(x),0,0,k(1),FUNC,REAL,false), [-pi/2, pi/2], 'trig'), M(2));
        c12 = trigcoeffs(chebfun(@(x) u_exact(l(2),yt(x),1,0,k(1),FUNC,REAL,false), [-pi/2, pi/2], 'trig'), M(2));
    end

    if(SUBDOMAIN.basis(3) == 'c')
        c03 = chebcoeffs(chebfun(@(x) u_exact(xt(x),l(3),0,0,k(1),FUNC,REAL,false), [-1,1]), M(3));
        c13 = chebcoeffs(chebfun(@(x) -u_exact(xt(x),l(3),0,1,k(1),FUNC,REAL,false), [-1,1]), M(3));
    elseif(SUBDOMAIN.basis(3) == 'f')
        c03 = trigcoeffs(chebfun(@(x) u_exact(xt(x),l(3),0,0,k(1),FUNC,REAL,false), [-pi/2, pi/2], 'trig'), M(3));
        c13 = trigcoeffs(chebfun(@(x) -u_exact(xt(x),l(3),0,1,k(1),FUNC,REAL,false), [-pi/2, pi/2], 'trig'), M(3));
    end

    if(SUBDOMAIN.basis(4) == 'c')
        c04 = chebcoeffs(chebfun(@(x) u_exact(xt(x),l(4),0,0,k(1),FUNC,REAL,false), [-1,1]), M(4));
        c14 = chebcoeffs(chebfun(@(x) u_exact(xt(x),l(4),0,1,k(1),FUNC,REAL,false), [-1,1]), M(4));
    elseif(SUBDOMAIN.basis(4) == 'f')
        c04 = trigcoeffs(chebfun(@(x) u_exact(xt(x),l(4),0,0,k(1),FUNC,REAL,false), [-pi/2, pi/2], 'trig'), M(4));
        c14 = trigcoeffs(chebfun(@(x) u_exact(xt(x),l(4),0,1,k(1),FUNC,REAL,false), [-pi/2, pi/2], 'trig'), M(4));
    end
    
    
    c = [c01, c11; c02, c12; c03, c13; c04, c14];
    if(rad>0)        
        c05 = trigcoeffs(chebfun(@(s) u_exact(rad, s, 0, 0, k, ...
                                        FUNC, REAL, true), [-pi,pi], 'trig'), Mr);
        c15 = trigcoeffs(chebfun(@(s) u_exact(rad, s, 1, 0, k, ...
                                        FUNC, REAL, true), [-pi,pi], 'trig'), Mr);
                                    
        c = [c; [c05, c15]];
    end
end