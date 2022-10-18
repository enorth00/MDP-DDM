%%%---------------------------------------------%%%
% Function name: rotate_functions
% 
% Purpose: Accept information computed from basis functions along side 1, 
%          return that info directly rotated to corresponding 
%          positions on another side.
% 
% Arguments:  x - data to be transferred. If x is a matrix, then each
%                 column is individually operated upon.
%      side_out - which side the data should be transferred to
%         sides - a sides grid to aid transfer. 
%                 (Note: sides(sides>0) is equivalent to
%                  Ngamma(Ngamma) as a logical matrix)
% 
% 
% 
% 
%%%---------------------------------------------%%%
function out = rotate_functions(x, side_out, sides)
    [N1,N2] = size(sides);
    husk = zeros(N1,N2);
    switch(side_out)
        case 1
            k=0;
        case 2
            k=2;
        case 3
            k=1;
        case 4
            k=3;
    end
    
    [N,M] = size(x);
    out = zeros(N,M);
    for i=1:M
        husk(sides>0) = x(:,i);
        B = conj(rot90(husk,k)) * (-1)^(i+1);
        out(:,i) = B(sides>0);
    end
end