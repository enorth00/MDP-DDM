function v = extend_xi_inhomogeneous(BLOCK, SUBDOMAIN)

    s = SUBDOMAIN.source;
    f = s.f; 
    fdx = s.fdx; 
    fdy = s.fdy; 
    fdx2 = s.fdx2; 
    fdy2 = s.fdy2;
    
    k = SUBDOMAIN.params.k(1); rad = BLOCK.radius;
    coords = BLOCK.grids.coords(:,1:2); sides = BLOCK.grids.coords(:,3);
%     e = SUBDOMAIN.edges; 
    e = BLOCK.edges; 
    LM = [e, rad];
    LL = @(s) LM(s);
    bounds = LL(sides)';
    
    [dist, X] = distance_calculator(coords, sides, bounds);
    
    % Four (five with rod) cases, one for each side
    v2 = zeros(length(sides), 1);
    v3 = zeros(length(sides), 1);
    v4 = zeros(length(sides), 1);
    
    idx1 = sides==1;
%     v2(idx1) = evaluate(@(x) f(e(1),x,k), X(idx1));
%     v3(idx1) = evaluate(@(x) -fdx(e(1),x,k), X(idx1));
%     v4(idx1) = evaluate(@(x) fdx2(e(1),x,k) - fdy2(e(1),x,k) - k^2*f(e(1),x,k), X(idx1));
    v2(idx1) = f(e(1),X(idx1),k);
    v3(idx1) = -fdx(e(1),X(idx1),k);
    v4(idx1) = fdx2(e(1),X(idx1),k) - fdy2(e(1),X(idx1),k) - k^2*f(e(1),X(idx1),k);
    
    idx2 = sides==2;
%     v2(idx2) = evaluate(@(x) f(e(2),x,k), X(idx2));
%     v3(idx2) = evaluate(@(x) fdx(e(2),x,k), X(idx2));
%     v4(idx2) = evaluate(@(x) fdx2(e(2),x,k) - fdy2(e(2),x,k) - k^2*f(e(2),x,k), X(idx2));
    v2(idx2) = f(e(2),X(idx2),k);
    v3(idx2) = fdx(e(2),X(idx2),k);
    v4(idx2) = fdx2(e(2),X(idx2),k) - fdy2(e(2),X(idx2),k) - k^2*f(e(2),X(idx2),k);

    idx3 = sides==3;
    v2(idx3) = f(X(idx3),e(3),k);
    v3(idx3) = -fdy(X(idx3),e(3),k);
    v4(idx3) = fdy2(X(idx3),e(3),k) - fdx2(X(idx3),e(3),k) - k^2*f(X(idx3),e(3),k);
    
    idx4 = sides==4;
    v2(idx4) = f(X(idx4),e(4),k);
    v3(idx4) = fdy(X(idx4),e(4),k);
    v4(idx4) = fdy2(X(idx4),e(4),k) - fdx2(X(idx4),e(4),k) - k^2*f(X(idx4),e(4),k);
    
    if(rad>0)
        idx5 = sides==5;
        v2(idx5) = s.rf(rad,X(idx5),k);
        v3(idx5) = s.rfdr(rad,X(idx5),k) - (1/rad)*v2(idx5);
        v4(idx5) = s.rfdr2(rad,X(idx5),k) + (2/rad^2 - k^2)*v2(idx5) - ...
                   (1/rad)*v3(idx5) - (1/rad^2)*s.rfdt2(rad,X(idx5),k);
    end
    
    v = (v2/2).*dist.^2 + (v3/6).*dist.^3 + (v4/24).*dist.^4;
end

function [dist, X] = distance_calculator(coords, sides, bounds)
    % Get the number of points for evaluation
    N = length(sides);

    % Allocate variables:
    % For each node, Xd is the coordinate used for calculating distance, 
    % and X is the coordinate used to evaluate the boundary function.
    Xd = zeros(N,1);
    X = zeros(N,1);

    % For sides 1 and 2, assign the x-coordinates for distance and
    % y-coordinates for evaluation
    idx12 = sides==1|sides==2;
    Xd(idx12) = coords(idx12, 1);
    X(idx12) = coords(idx12, 2);
    
    % For sides 3 and 4, assign the y-coordinates for distance and
    % x-coordinates for evaluation
    idx34 = sides==3|sides==4;
    Xd(idx34) = coords(idx34, 2);
    X(idx34) = coords(idx34, 1);
    
    idx5 = sides==5;
    [X(idx5), Xd(idx5)] = cart2pol(coords(idx5, 1), ...
                                           coords(idx5,2));
    
    dist = Xd - bounds;
    
    % For sides 1 and 3, invert the sign to reflect a positive outward
    % distance and a negative inward distance.
    dist(sides==1 | sides==3) = ...
        dist(sides==1 | sides==3) * (-1);
end