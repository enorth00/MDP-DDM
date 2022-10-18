% Returns a whole column of \gamma nodes that have been extended

%%% coords is the vector of (x,y) values in \gamma
%%% sides is the vector that indicates which side each element of coords
% belongs to.
%%% 

function v = extend_square_basisfuncs(xi0, xi1, xi0d2, xi1d2, xi0d4, r, BLOCK, SUBDOMAIN, k, P)

    coords = BLOCK.grids.coords(:,1:2);
    sides = BLOCK.grids.coords(:,3);
    
%     k = SUBDOMAIN.params.k;
    M = SUBDOMAIN.params.M(1:4);
    Mr = SUBDOMAIN.params.M(5);
    rad = BLOCK.radius;
    L = BLOCK.edges;
    
    a = L(1); b = L(2); c = L(3); d = L(4);
    
    xt = @(x) (x - (a+b)/2) * (2/(b-a));
    yt = @(y) (y - (c+d)/2) * (2/(d-c));
    
%     x_bar = coords(1); y_bar = coords(2);
    x = xt(coords(:,1));
    y = yt(coords(:,2));
    
    % convert "sides" to the actual coordinate
    LM = [L, rad];
    LL = @(s) LM(s);
    bounds = LL(sides)';
    
    [dist,X] = distance_calculator([x,y], sides, bounds);
%     distv = repmat(dist,1,M); 
    distr = repmat(dist(sides==5),1,Mr);
    Z = zeros(length(coords(sides==5)), Mr);
    
    if(P==0)
        v0 = xi0(X,:);
        v2 = -xi0d2(X,:) - k^2*v0;
        v4 = xi0d4(X,:) + 2*k^2*xi0d2(X,:) + k^4*v0;
        
        V = v0 + v2.*(dist.^2)/2 + v4.*(dist.^4)/24; 
        if(Mr>0)
            
            idx = sides==5;
            v0r = r.d0.xi0(X(idx),:);
            v1r = Z;
            v2r = -k^2*v0r - (1/rad^2)*r.d2.xi0(X(idx),:) - (1/rad)*v1r;
            v3r = (-k^2)*v1r + (2/rad^3)*r.d2.xi0(X(idx),:) + ...
                  (1/rad^2)*v1r - (1/rad)*v2r;
            v4r = (-2/rad^3)*v1r + (2/rad^2 - k^2)*v2r - (1/rad)*v3r + ...
                  (k^2/rad^2 - 6/rad^4)*r.d2.xi0(X(idx),:) + ...
                  (1/rad^4)*r.d4.xi0(X(idx),:);

            Vr = v0r + v1r.*distr + v2r.*(distr.^2/2) + ...
                              v3r.*(distr.^3/6) + v4r.*(distr.^4/24);
            
%             idx = sides==5;
%             v0r = r.d0.xi0(X(idx),:);
%             v2r = -k^2*v0r - (1/rad^2)*r.d2.xi0(X(idx)); 
%             v3r = -(2/rad^3)*r.d2.xi0(X(idx)) - (1/rad)*v2r;
%             v4r = (k^4 - 3*k^2/rad^2)*v0r + ...
%                   (2*k^2/rad^2 - 11/rad^4)*r.d2.xi0(X(idx)) + ...
%                   (1/rad^4)*r.d4.xi0(X(idx));
% 
%             Vr = v0r + v2r.*(distr.^2/2) + v3r.*(distr.^3/6) + ...
%                     v4r.*(distr.^4/24); 
%             idx = sides==5;
%             v0r = r.d0.xi0(X(idx));
% %             v1r = r.d0.xi1(X(idx));
%             v1r = Z;
%             v2r = -k^2*v0r - (1/rad^2)*r.d2.xi0(X(idx)) - (1/rad)*v1r;
%             v3r = (-k^2)*v1r - (2/rad^3)*r.d2.xi0(X(idx)) + ...
%                     (1/rad^2)*v1r + (1/rad)*v2r;
%             v4r = (-2/rad^3)*v1r + (2/rad^2 - k^2)*v2r - (1/rad)*v3r + ...
%               (k^2/rad^2 - 6/rad^4)*r.d2.xi0(X(idx)) + ...
%               (1/rad^4)*r.d4.xi0(X(idx));
% 
%             Vr = v0r + v1r.*distr + v2r.*(distr.^2/2) + ...
%                               v3r.*(distr.^3/6) + v4r.*(distr.^4/24);
        end  
    else
        
        v1 = xi1(X,:);
        v3 = -xi1d2(X,:) - k^2*v1;

        V = v1.*dist + v3.*(dist.^3)/6;

        
        if(Mr>0)
            
            idx = sides==5;
            v0r = Z;
            v1r = r.d0.xi1(X(idx),:);
            v2r = -k^2*v0r - (1/rad)*v1r;
            v3r = (-k^2)*v1r - (1/rad^2)*r.d2.xi1(X(idx),:) + ...
                  (1/rad^2)*v1r - (1/rad)*v2r;
            v4r = (-2/rad^3)*v1r + (2/rad^2 - k^2)*v2r - (1/rad)*v3r + ...
                  (5/rad^3)*r.d2.xi1(X(idx),:);

            Vr = v0r + v1r.*distr + v2r.*(distr.^2/2) + ...
                              v3r.*(distr.^3/6) + v4r.*(distr.^4/24);
            
%             idx = sides==5;
% %             v0r = r.d0.xi0(X(idx));
%             v0r = Z;
%             v1r = r.d0.xi1(X(idx));
%             v2r = -k^2*v0r - (1/rad)*v1r;
%             v3r = (-k^2 + 1/rad^2)*v1r - (1/rad^2)*r.d2.xi1(X(idx)) - (1/rad)*v2r;
%             v4r = (-2/rad^3)*v1r + (2/rad^2 - k^2)*v2r - (1/rad)*v3r + ...
%                     (5/rad^3)*r.d2.xi1(X(idx));
% 
%             Vr = v0r + v1r.*distr(idx) + v2r.*(distr(idx).^2/2) + ...
%                               v3r.*(distr(idx).^3/6) + v4r.*(distr(idx).^4/24);
        end
    end
    
    NN = length(sides);
    v = zeros(NN, sum(M)+Mr);
    
  %%% Distribute the node evaluations to the appropriate columns of the
  %%% answer matrix.
    v(sides==1,1:M) = V(sides==1,:);
    v(sides==2,M+1:2*M) = V(sides==2,:);
    v(sides==3,2*M+1:3*M) = V(sides==3,:);
    v(sides==4,3*M+1:4*M) = V(sides==4,:);
    if(BLOCK.radius>0)
        v(sides==5,4*M+1:4*M+Mr) = Vr;
    end
    
end


function [dist, X] = distance_calculator(coords, sides, bounds)
    N = length(sides);
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