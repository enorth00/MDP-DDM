% sides - vector obtained with gamma_trace
% xi structures - each contains xi0_1, xi1_1, xi0_2, etc...

function xi_gamma = extend_xi(d0,d2,d4,r, BLOCK, SUBDOMAIN, F)
    
    % Convenience defs for several variables
    coords = BLOCK.grids.coords(:,1:2);
    sides = BLOCK.grids.coords(:,3);
    k = SUBDOMAIN.params.k(1);
    rad = BLOCK.radius;

    % convert "sides" to the actual coordinate
    LM = [BLOCK.edges, rad];
    LL = @(s) LM(s);
    bounds = LL(sides)';

    % dist contains the coordinate used for that side's distance
    % calculation, and X contains the coordinate used in evaluation of the
    % boundary function.
    [dist,X] = distance_calculator(coords, sides, bounds);

    % Compute each member of the Taylor-style extension. We do this for
    % EVERY node at EVERY function because vectorizing is faster than doing
    % it node-by-node.
    v0 = zeros(length(X),4); v1 = v0; v2 = v0; v3 = v0; v4 = v0;
    for i=1:4
        v0(:,i) = d0.xi0{i}(X);
        v1(:,i) = d0.xi1{i}(X);
        v2(:,i) = -d2.xi0{i}(X) - k^2*v0(:,i);
        v3(:,i) = -d2.xi1{i}(X) - k^2*v1(:,i);
        v4(:,i) = d4.xi0{i}(X) + 2*k^2*d2.xi0{i}(X) + k^4*v0(:,i);
    end
    
%     v0 = d0.xi0(X,:);
%     v1 = d0.xi1(X,:);
%     v2 = -d2.xi0(X,:) - k^2*v0;
%     v3 = -d2.xi1(X,:) - k^2*v1;
%     v4 = d4.xi0(X,:) + 2*k^2*d2.xi0(X,:) + k^4*v0;

    XI_gamma_h = v0 + v1.*dist + v2.*(dist.^2/2) + ...
                 v3.*(dist.^3/6) + v4.*(dist.^4/24);

    % If necessary, also do the polar ones for the rod.
    if(rad>0)
        idx = sides==5;
        v0r = r.d0.xi0{1}(X(idx));
        v1r = r.d0.xi1{1}(X(idx));
        v2r = -k^2*v0r - (1/rad^2)*r.d2.xi0{1}(X(idx)) - (1/rad)*v1r;
        v3r = (-k^2)*v1r + (2/rad^3)*r.d2.xi0{1}(X(idx)) - ...
              (1/rad^2)*r.d2.xi1{1}(X(idx)) + (1/rad^2)*v1r - (1/rad)*v2r;
        v4r = (-2/rad^3)*v1r + (2/rad^2 - k^2)*v2r - (1/rad)*v3r + ...
              (k^2/rad^2 - 6/rad^4)*r.d2.xi0{1}(X(idx)) + ...
              (5/rad^3)*r.d2.xi1{1}(X(idx)) + (1/rad^4)*r.d4.xi0{1}(X(idx));
                
        XI_gamma_R = v0r + v1r.*dist(idx) + v2r.*(dist(idx).^2/2) + ...
                          v3r.*(dist(idx).^3/6) + v4r.*(dist(idx).^4/24);

    end
    
    % Split information into appropriate elements of a final solution
    % vector.
    xi_gamma_h = zeros(length(sides),1);
    xi_gamma_h(sides==1) = XI_gamma_h(sides==1,1);
    xi_gamma_h(sides==2) = XI_gamma_h(sides==2,2);
    xi_gamma_h(sides==3) = XI_gamma_h(sides==3,3);
    xi_gamma_h(sides==4) = XI_gamma_h(sides==4,4);
    if(BLOCK.radius>0)
        xi_gamma_h(sides==5) = XI_gamma_R;
    end

    % Add on inhomogeneous source information
    xi_gamma = xi_gamma_h + F;
    
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