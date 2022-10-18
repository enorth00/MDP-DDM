% Returns a whole column of \gamma nodes that have been extended

%%% coords is the vector of (x,y) values in \gamma
%%% sides is the vector that indicates which side each element of coords
% belongs to.
%%% 

function Vr = extend_rod_basisfuncs(r, BLOCK, SUBDOMAIN, P)

    coords1 = BLOCK.grids.coords(:,1:2);
    sides = BLOCK.grids.coords(:,3);

    coords = coords1(sides==5,:);

    k = SUBDOMAIN.params.k(2);
    [~,Mr] = size(r.d0.xi0);
    rad = BLOCK.radius;
    L = BLOCK.edges;
    
    a = L(1); b = L(2); c = L(3); d = L(4);
    
    xt = @(x) (x - (a+b)/2) * (2/(b-a));
    yt = @(y) (y - (c+d)/2) * (2/(d-c));
    
    x = xt(coords(:,1));
    y = yt(coords(:,2));
    
    [dist,X] = distance_calculator([x,y], rad);
    distr = repmat(dist,1,Mr);
    Z = zeros(length(coords), Mr);
    
    if(P==0)
        v0r = r.d0.xi0{1}(X,:);
        v1r = Z;
        v2r = -k^2*v0r - (1/rad^2)*r.d2.xi0{1}(X,:) - (1/rad)*v1r;
        v3r = (-k^2)*v1r + (2/rad^3)*r.d2.xi0{1}(X,:) + ...
              (1/rad^2)*v1r - (1/rad)*v2r;
        v4r = (-2/rad^3)*v1r + (2/rad^2 - k^2)*v2r - (1/rad)*v3r + ...
              (k^2/rad^2 - 6/rad^4)*r.d2.xi0{1}(X,:) + ...
              (1/rad^4)*r.d4.xi0{1}(X,:);

        Vr = v0r + v1r.*distr + v2r.*(distr.^2/2) + ...
                          v3r.*(distr.^3/6) + v4r.*(distr.^4/24);
    else
        v0r = Z;
        v1r = r.d0.xi1{1}(X,:);
        v2r = -k^2*v0r - (1/rad)*v1r;
        v3r = (-k^2)*v1r - (1/rad^2)*r.d2.xi1{1}(X,:) + ...
              (1/rad^2)*v1r - (1/rad)*v2r;
        v4r = (-2/rad^3)*v1r + (2/rad^2 - k^2)*v2r - (1/rad)*v3r + ...
              (5/rad^3)*r.d2.xi1{1}(X,:);

        Vr = v0r + v1r.*distr + v2r.*(distr.^2/2) + ...
                          v3r.*(distr.^3/6) + v4r.*(distr.^4/24);
    end
    
end


function [dist, X] = distance_calculator(coords, rad)
    bounds = rad*ones(length(coords(:,1)),1);
    [X, Xd] = cart2pol(coords(:, 1), coords(:,2));
    
    dist = Xd - bounds;
end