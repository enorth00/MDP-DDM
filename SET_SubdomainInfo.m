function SUBDOMAIN = SET_SubdomainInfo(SUBDOMAIN, BC, ADJ, OPTIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   SET_SubdomainInfo
%
%   Inputs:     
%               
%           
%   Purpose:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SUBDOMAIN.boundary = set_boundary(SUBDOMAIN, BC, ADJ, OPTIONS);
    SUBDOMAIN.source = set_source(SUBDOMAIN, OPTIONS);
end

function s = set_source(SUBDOMAIN, OPTIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   set_source
%
%   Inputs:     
%               
%           
%   Purpose:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FUNC = OPTIONS.FUNC;
    REAL = OPTIONS.REAL;
    
%     k = SUBDOMAIN.params.k;
    
    a = (SUBDOMAIN.edges(2) + SUBDOMAIN.edges(1))/2;
    b = (SUBDOMAIN.edges(4) + SUBDOMAIN.edges(3))/2;
 
    s.f = @(x,y,k) source(x,y,0,0,k,FUNC,REAL,false,a,b);
    s.fdx = @(x,y,k) source(x,y,1,0,k,FUNC,REAL,false,a,b);
    s.fdy = @(x,y,k) source(x,y,0,1,k,FUNC,REAL,false,a,b);
    s.fdx2 = @(x,y,k) source(x,y,2,0,k,FUNC,REAL,false,a,b);
    s.fdy2 = @(x,y,k) source(x,y,0,2,k,FUNC,REAL,false,a,b);
        
    s.rf = @(r,t,k) source(r,t,0,0,k,FUNC,REAL,true,a,b);
    s.rfdr = @(r,t,k) source(r,t,1,0,k,FUNC,REAL,true,a,b);
    s.rfdt = @(r,t,k) source(r,t,0,1,k,FUNC,REAL,true,a,b);
    s.rfdr2 = @(r,t,k) source(r,t,2,0,k,FUNC,REAL,true,a,b);
    s.rfdt2 = @(r,t,k) source(r,t,0,2,k,FUNC,REAL,true,a,b);

    s.FUNC = FUNC;
end

function boundary = set_boundary(SUBDOMAIN, BC, adj, OPTIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   set_boundary
%
%   Inputs:     SUBDOMAIN - 
%                      BC -
%                     adj - 
%                 OPTIONS - 
%           
%   Purpose:    Sets the function handles of each boundary condition along
%               the corresponding edge in the subdomain. Handles are stored
%               as chebfuns (either cheb or trig).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u = @(x,y,dx,dy) u_exact(x,y,dx,dy,SUBDOMAIN.params.k(1),...
                             OPTIONS.FUNC,OPTIONS.REAL,false);
    boundary = struct('fun', cell(1,4), ...
                   'alpha', cell(1,4), ...
                   'beta', cell(1,4)...
                   );
    a = BC(1,:); b = BC(2,:);
    for i=1:4
        if(adj(i) ~= 0)
            a(i) = 0;
            b(i) = 0;
        end
    end
    
    for i=1:4
        if(SUBDOMAIN.basis(i) == 'c')
            boundary(i) = chebyshev_basis(SUBDOMAIN, [a(i); b(i)], u, i);
        elseif(SUBDOMAIN.basis(i) == 'f')
            boundary(i) = fourier_basis(SUBDOMAIN, [a(i); b(i)], u, i, OPTIONS.FOURIER_AUX);
        end
    end
end

function bound = chebyshev_basis(SUBDOMAIN,BCs,u,side)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   chebyshev_basis
%
%   Inputs:     
%               
%           
%   Purpose:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = BCs(1); b = BCs(2);
    l = SUBDOMAIN.edges;
    
    switch(side)
        case 1
            bound.fun = chebfun(@(x) a*u(l(1),x,0,0) - b*u(l(1),x,1,0), [l(3),l(4)]);
        case 2
            bound.fun = chebfun(@(x) a*u(l(2),x,0,0) + b*u(l(2),x,1,0), [l(3),l(4)]);
        case 3
            bound.fun = chebfun(@(x) a*u(x,l(3),0,0) - b*u(x,l(3),0,1), [l(1),l(2)]);
        case 4
            bound.fun = chebfun(@(x) a*u(x,l(4),0,0) + b*u(x,l(4),0,1), [l(1),l(2)]);
    end
    bound.alpha = a; bound.beta = b;
end
function bound = fourier_basis(SUBDOMAIN,BCs,u,side,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   fourier_basis
%
%   Inputs:     
%               
%           
%   Purpose:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = BCs(1); b = BCs(2); 
    extension_parameter = P(2);
    le = SUBDOMAIN.edges; 
    p = (extension_parameter - 1);
    l = [le(1)-p, le(2)+p, le(3)-p, le(4)+p];
    
    switch(side)
        case 1
%             X = [le(3:4), l(3:4)];
%             ur = chebfun(@(x) a*u(le(1),x,0,0) - b*u(le(1),x,1,0), X(1:2));
%             U = periodic_extension(ur, X);
%             bound.fun = chebfun(@(x) U(x), [l(3),l(4)], 'trig');
            bound.fun = chebfun(@(x) a*u(le(1),x,0,0) - b*u(le(1),x,1,0), [l(3),l(4)], 'trig');
        case 2
%             X = [le(3:4), l(3:4)];
%             ur = chebfun(@(x) a*u(le(2),x,0,0) + b*u(le(2),x,1,0), X(1:2));
%             U = periodic_extension(ur, X);
%             bound.fun = chebfun(@(x) U(x), [l(3),l(4)], 'trig');
            bound.fun = chebfun(@(x) a*u(le(2),x,0,0) + b*u(le(2),x,1,0), [l(3),l(4)], 'trig');
        case 3
%             X = [le(1:2), l(1:2)];
%             ur = chebfun(@(x) a*u(x,le(3),0,0) - b*u(x,le(3),0,1), X(1:2));
%             U = periodic_extension(ur, X);
%             bound.fun = chebfun(@(x) U(x), [l(1),l(2)], 'trig');
            bound.fun = chebfun(@(x) a*u(x,le(3),0,0) - b*u(x,le(3),0,1), [l(1),l(2)], 'trig');
        case 4
%             X = [le(1:2), l(1:2)];
%             ur = chebfun(@(x) a*u(x,le(4),0,0) + b*u(x,le(4),0,1), X(1:2));
%             U = periodic_extension(ur, X);
%             bound.fun = chebfun(@(x) U(x), [l(1),l(2)], 'trig');
            bound.fun = chebfun(@(x) a*u(x,le(4),0,0) + b*u(x,le(4),0,1), [l(1),l(2)], 'trig');
    end
    bound.alpha = a; bound.beta = b;
end