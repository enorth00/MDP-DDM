
% Debugging purposes.
% 
% Provides exact boundary function handles instead of reconstructing the
% handles from the coefficients.

function [d0, d2, d4, r] = DEBUG_get_exact_boundary_functions(BLOCK, SUBDOMAIN, OPTIONS)
    e = SUBDOMAIN.edges; e1 = BLOCK.edges;
    a = e(1); b = e(2); c = e(3); d = e(4);
    a1 = e1(1); c1 = e1(3);
    k = SUBDOMAIN.params.k; FUNC = OPTIONS.FUNC;
    rad = BLOCK.radius;
    
    tmp = chebpoly(1:4);
    d0.xi0 = tmp;
    d0.xi1 = tmp;
    d2.xi0 = tmp;
    d2.xi1 = tmp;
    d4.xi0 = tmp;

    xt = @(x) x + (a - a1);
    yt = @(y) y + (c - c1);
    
    d0.xi0(:,1) = chebfun(@(x) u_exact(a,yt(x),0,0,k,FUNC,OPTIONS.REAL,false), [-1,1]);
    d0.xi1(:,1) = chebfun(@(x) -u_exact(a,yt(x),1,0,k,FUNC,OPTIONS.REAL,false), [-1,1]);

    d0.xi0(:,2) = chebfun(@(x) u_exact(b,yt(x),0,0,k,FUNC,OPTIONS.REAL,false), [-1,1]);
    d0.xi1(:,2) = chebfun(@(x) u_exact(b,yt(x),1,0,k,FUNC,OPTIONS.REAL,false), [-1,1]);
    
    d0.xi0(:,3) = chebfun(@(x) u_exact(xt(x),c,0,0,k,FUNC,OPTIONS.REAL,false), [-1,1]);
    d0.xi1(:,3) = chebfun(@(x) -u_exact(xt(x),c,0,1,k,FUNC,OPTIONS.REAL,false), [-1,1]);
    
    d0.xi0(:,4) = chebfun(@(x) u_exact(xt(x),d,0,0,k,FUNC,OPTIONS.REAL,false), [-1,1]);
    d0.xi1(:,4) = chebfun(@(x) u_exact(xt(x),d,0,1,k,FUNC,OPTIONS.REAL,false), [-1,1]);
    
    for i=1:4
        d2.xi0(:,i) = diff(d0.xi0(:,i), 2);
        d2.xi1(:,i) = diff(d0.xi1(:,i), 2);
        
        d4.xi0(:,i) = diff(d2.xi0(:,i), 2);
    end
    
    
    if(rad>0)
        r.d0.xi0 = chebfun(@(s) u_exact(rad,s,0,0,k,FUNC,OPTIONS.REAL,true), [-pi,pi], 'trig');
        r.d0.xi1 = chebfun(@(s) u_exact(rad,s,1,0,k,FUNC,OPTIONS.REAL,true), [-pi,pi], 'trig');
        
        r.d2.xi0 = diff(r.d0.xi0, 2);
        r.d2.xi1 = diff(r.d0.xi1, 2);
        
        r.d4.xi0 = diff(r.d2.xi0, 2);
    else
        r = 0;
    end
end