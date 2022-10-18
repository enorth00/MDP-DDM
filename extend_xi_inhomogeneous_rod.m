function v = extend_xi_inhomogeneous_rod(BLOCK, SUBDOMAIN)

    s = SUBDOMAIN.source;
    
    k = SUBDOMAIN.params.k(2); 
    rad = BLOCK.radius;
    
    sides = BLOCK.grids.coords(:,3); 
    coords = BLOCK.grids.coords(sides==5,1:2); 
    
    [dist, X] = distance_calculator(coords, rad);

    v2 = s.rf(rad,X,k);
    v3 = s.rfdr(rad,X,k) - (1/rad)*v2;
    v4 = s.rfdr2(rad,X,k) + (2/rad^2 - k^2)*v2 - (1/rad)*v3 - ...
            (1/rad^2)*s.rfdt2(rad,X,k);
    
    v = (v2/2).*dist.^2 + (v3/6).*dist.^3 + (v4/24).*dist.^4;
end

function [dist, X] = distance_calculator(coords, rad)
    bounds = rad*ones(length(coords(:,1)),1);
    [X, Xd] = cart2pol(coords(:, 1), coords(:,2));
    
    dist = Xd - bounds;
end