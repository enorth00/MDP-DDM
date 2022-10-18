%-------------------------------------------------------------------------
%
%   Function: source
%
%   Inputs:     x - x value for evaluation
%               y - y value for evaluation
%              dx - number of x derivatives
%              dy - number of y derivatives
%               r - which test source to use
%                      0 - custom
%                      1 ... 8 - preset, see 'test_source.m'
%               k - wavenumber
%
%   Output:     cg - indicated function evaluated at (x,y)
%
%   Description:    This function provides a central location from which
%                   from which to obtain the exact source, no matter
%                   which case is being tested.
%   
%-------------------------------------------------------------------------

function cg = source(x,y,dx,dy,k,r,REAL,POL,a,b)
    
    if(exist('POL', 'var') == 0)
        POL = false;
    end
    if (r==0)
        % Placeholder for a future possibility
        cg = 0;
    elseif(POL)
        rad = x; t = y;
        [x,y] = pol2cart(t,rad);
        cg = 0;
        if(dx==0 && dy==0)
            cg = test_source(x+a,y+b,0,0,k,r,false);
        elseif(dx==1 && dy==0)
            % first "r" derivative
            cg = test_source(x+a,y+b,1,0,k,r,false) .* cos(t) + ...
                 test_source(x+a,y+b,0,1,k,r,false) .* sin(t);
        elseif(dx==2 && dy==0)
            % second "r" derivative
            cg = test_source(x+a,y+b,2,0,k,r,false) .* (cos(t)).^2 + ...
                 test_source(x+a,y+b,0,2,k,r,false) .* (sin(t)).^2 + ...
                 2*test_source(x+a,y+b,1,1,k,r,false) .* cos(t) .* sin(t);
        elseif(dx==0 && dy==1)
            cg = test_source(x+a,y+b,1,0,k,r,false) .* (-rad.*sin(t)) + ...
                 test_source(x+a,y+b,0,1,k,r,false) .* (rad.*cos(t));
        elseif(dx==0 && dy==2)
            cg = test_source(x+a,y+b,2,0,k,r,false) .* (-rad.*sin(t)).^2 + ...
                 test_source(x+a,y+b,0,2,k,r,false) .* (rad.*cos(t)).^2 + ...
                 test_source(x+a,y+b,1,1,k,r,false) .* 2.*(-rad.^2.*cos(t).*sin(t)) + ...
                 test_source(x+a,y+b,1,0,k,r,false) .* (-rad.*cos(t)) + ...
                 test_source(x+a,y+b,0,1,k,r,false) .* (-rad.*sin(t));
        end
    else
        cg = test_source(x+a,y+b,dx,dy,k,r,POL);
    end
    if(REAL)
        cg = real(cg);
    end
end

function cg = test_source(x,y,dx,dy,k,r,POL)
    k = k(1);
    switch r
        case 1
            % homogeneous case, no source
            cg = 0;
            
        case 2
            cg = 0;
            
        case 3
            cg = 0;
           
            
            
            
            
        case 4               
            if(x.^2 + y.^2 >= 1)
                cg = 0;
            elseif(dx==0 && dy==0)
                cg = (k .^ 2 .* (x .^ 2 + y .^ 2 - 1) .^ 4 + 4 .* x .^ 4 + (8 .* y .^ 2 + 4) .* x .^ 2 + 4 .* y .^ 4 + 4 .* y .^ 2 - 4) .* exp(0.1e1 ./ (x .^ 2 + y .^ 2 - 1)) ./ ((x .^ 2 + y .^ 2 - 1) .^ 4);
            elseif(dx==1 && dy==0)
                cg = -0.2e1 .* x .* (k .^ 2 .* x .^ 8 + (4 .* k .^ 2 .* y .^ 2 - 4 .* k .^ 2 + 8) .* x .^ 6 + (6 .* k .^ 2 .* y .^ 4 + (-12 .* k .^ 2 + 24) .* y .^ 2 + 6 .* k .^ 2 + 16) .* x .^ 4 + (4 .* k .^ 2 .* y .^ 6 + (-12 .* k .^ 2 + 24) .* y .^ 4 + (12 .* k .^ 2 + 32) .* y .^ 2 - 4 .* k .^ 2 - 28) .* x .^ 2 + k .^ 2 .* y .^ 8 + (-4 .* k .^ 2 + 8) .* y .^ 6 + (6 .* k .^ 2 + 16) .* y .^ 4 + (-4 .* k .^ 2 - 28) .* y .^ 2 + k .^ 2 + 8) .* exp(0.1e1 ./ (x .^ 2 + y .^ 2 - 1)) ./ ((x .^ 2 + y .^ 2 - 1) .^ 6);
            elseif(dx==0 && dy==1)
                cg = -0.2e1 .* y .* (k .^ 2 .* y .^ 8 + (4 .* k .^ 2 .* x .^ 2 - 4 .* k .^ 2 + 8) .* y .^ 6 + (6 .* k .^ 2 .* x .^ 4 + (-12 .* k .^ 2 + 24) .* x .^ 2 + 6 .* k .^ 2 + 16) .* y .^ 4 + (4 .* k .^ 2 .* x .^ 6 + (-12 .* k .^ 2 + 24) .* x .^ 4 + (12 .* k .^ 2 + 32) .* x .^ 2 - 4 .* k .^ 2 - 28) .* y .^ 2 + k .^ 2 .* x .^ 8 + (-4 .* k .^ 2 + 8) .* x .^ 6 + (6 .* k .^ 2 + 16) .* x .^ 4 + (-4 .* k .^ 2 - 28) .* x .^ 2 + k .^ 2 + 8) .* exp(0.1e1 ./ (x .^ 2 + y .^ 2 - 1)) ./ ((x .^ 2 + y .^ 2 - 1) .^ 6);
            elseif(dx==2 && dy==0)
                cg = 0.6e1 .* (k .^ 2 .* x .^ 12 + (0.40e2 ./ 0.3e1 + 0.14e2 ./ 0.3e1 .* k .^ 2 .* y .^ 2 - 0.4e1 .* k .^ 2) .* x .^ 10 + (0.25e2 ./ 0.3e1 .* k .^ 2 .* y .^ 4 + (-0.14e2 .* k .^ 2 + 0.152e3 ./ 0.3e1) .* y .^ 2 + 0.17e2 ./ 0.3e1 .* k .^ 2 + 0.48e2) .* x .^ 8 + (0.20e2 ./ 0.3e1 .* k .^ 2 .* y .^ 6 + (-0.16e2 .* k .^ 2 + 0.208e3 ./ 0.3e1) .* y .^ 4 + (0.12e2 .* k .^ 2 + 0.144e3) .* y .^ 2 - 0.8e1 ./ 0.3e1 .* k .^ 2 - 0.308e3 ./ 0.3e1) .* x .^ 6 + (0.5e1 ./ 0.3e1 .* k .^ 2 .* y .^ 8 + (-0.4e1 .* k .^ 2 + 0.112e3 ./ 0.3e1) .* y .^ 6 + (0.2e1 .* k .^ 2 + 0.144e3) .* y .^ 4 + (-0.188e3 + 0.4e1 ./ 0.3e1 .* k .^ 2) .* y .^ 2 - k .^ 2 + 0.40e2) .* x .^ 4 + (-0.2e1 ./ 0.3e1 .* k .^ 2 .* y .^ 10 + (0.4e1 .* k .^ 2 + 0.8e1 ./ 0.3e1) .* y .^ 8 + (-0.28e2 ./ 0.3e1 .* k .^ 2 + 0.48e2) .* y .^ 6 + (-0.68e2 + 0.32e2 ./ 0.3e1 .* k .^ 2) .* y .^ 4 + (0.40e2 ./ 0.3e1 - 0.6e1 .* k .^ 2) .* y .^ 2 + 0.4e1 ./ 0.3e1 .* k .^ 2 + 0.20e2 ./ 0.3e1) .* x .^ 2 - (y - 0.1e1) .^ 2 .* (k .^ 2 .* y .^ 8 + (-0.4e1 .* k .^ 2 + 0.8e1) .* y .^ 6 + (0.6e1 .* k .^ 2 + 0.16e2) .* y .^ 4 + (-0.4e1 .* k .^ 2 - 0.28e2) .* y .^ 2 + k .^ 2 + 0.8e1) .* (y + 0.1e1) .^ 2 ./ 0.3e1) .* exp(0.1e1 ./ (x .^ 2 + y .^ 2 - 0.1e1)) ./ (x .^ 2 + y .^ 2 - 0.1e1) .^ 8;
            elseif(dx==0 && dy==2)
                cg = -0.2e1 .* (-3 .* k .^ 2 .* y .^ 12 + (-14 .* k .^ 2 .* x .^ 2 + 12 .* k .^ 2 - 40) .* y .^ 10 + (-25 .* k .^ 2 .* x .^ 4 + (42 .* k .^ 2 - 152) .* x .^ 2 - 17 .* k .^ 2 - 144) .* y .^ 8 + (-20 .* k .^ 2 .* x .^ 6 + (48 .* k .^ 2 - 208) .* x .^ 4 + (-36 .* k .^ 2 - 432) .* x .^ 2 + 8 .* k .^ 2 + 308) .* y .^ 6 + (-5 .* k .^ 2 .* x .^ 8 + (12 .* k .^ 2 - 112) .* x .^ 6 + (-6 .* k .^ 2 - 432) .* x .^ 4 + (-4 .* k .^ 2 + 564) .* x .^ 2 + 3 .* k .^ 2 - 120) .* y .^ 4 + (2 .* k .^ 2 .* x .^ 10 + (-12 .* k .^ 2 - 8) .* x .^ 8 + (28 .* k .^ 2 - 144) .* x .^ 6 + (-32 .* k .^ 2 + 204) .* x .^ 4 + (18 .* k .^ 2 - 40) .* x .^ 2 - 4 .* k .^ 2 - 20) .* y .^ 2 + (k .^ 2 .* x .^ 8 + (-4 .* k .^ 2 + 8) .* x .^ 6 + (6 .* k .^ 2 + 16) .* x .^ 4 + (-4 .* k .^ 2 - 28) .* x .^ 2 + k .^ 2 + 8) .* (x - 1) .^ 2 .* (x + 1) .^ 2) .* exp(0.1e1 ./ (x .^ 2 + y .^ 2 - 1)) ./ ((x .^ 2 + y .^ 2 - 1) .^ 8);
            elseif(dx==1 && dy==1)
                cg = 0.8e1 .* x .* y .* (k .^ 2 .* x .^ 10 + (-0.9e1 ./ 0.2e1 .* k .^ 2 + 0.12e2 + 0.5e1 .* k .^ 2 .* y .^ 2) .* x .^ 8 + (0.10e2 .* k .^ 2 .* y .^ 4 + (-0.18e2 .* k .^ 2 + 0.48e2) .* y .^ 2 + 0.8e1 .* k .^ 2 + 0.36e2) .* x .^ 6 + (0.10e2 .* k .^ 2 .* y .^ 6 + (-0.27e2 .* k .^ 2 + 0.72e2) .* y .^ 4 + (0.24e2 .* k .^ 2 + 0.108e3) .* y .^ 2 - 0.7e1 .* k .^ 2 - 0.90e2) .* x .^ 4 + (0.5e1 .* k .^ 2 .* y .^ 8 + (-0.18e2 .* k .^ 2 + 0.48e2) .* y .^ 6 + (0.24e2 .* k .^ 2 + 0.108e3) .* y .^ 4 + (-0.14e2 .* k .^ 2 - 0.180e3) .* y .^ 2 + 0.3e1 .* k .^ 2 + 0.50e2) .* x .^ 2 + k .^ 2 .* y .^ 10 + (-0.9e1 ./ 0.2e1 .* k .^ 2 + 0.12e2) .* y .^ 8 + (0.8e1 .* k .^ 2 + 0.36e2) .* y .^ 6 + (-0.7e1 .* k .^ 2 - 0.90e2) .* y .^ 4 + (0.3e1 .* k .^ 2 + 0.50e2) .* y .^ 2 - 0.6e1 - k .^ 2 ./ 0.2e1) .* exp(0.1e1 ./ (x .^ 2 + y .^ 2 - 0.1e1)) ./ (x .^ 2 + y .^ 2 - 0.1e1) .^ 8;
            end
            
            
            
            
            
            
        case 5
            cg=0;
            
        case 6
            cg = 0;
            if(dx==0 && dy==0)
                cg = -0.3e1 .* exp(-x .^ 2) .* cos(y) + 0.4e1 .* x .^ 2 .* exp(-x .^ 2) .* cos(y) + k .^ 2 .* exp(-x .^ 2) .* cos(y);
            elseif(dx==1 && dy==0)
                cg = 0.14e2 .* x .* exp(-x .^ 2) .* cos(y) - 0.8e1 .* x .^ 3 .* exp(-x .^ 2) .* cos(y) - 0.2e1 .* k .^ 2 .* x .* exp(-x .^ 2) .* cos(y);
            elseif(dx==0 && dy==1)
                cg = 0.3e1 .* exp(-x .^ 2) .* sin(y) - 0.4e1 .* x .^ 2 .* exp(-x .^ 2) .* sin(y) - k .^ 2 .* exp(-x .^ 2) .* sin(y);
            elseif(dx==2 && dy==0)
                cg = 0.14e2 .* exp(-x .^ 2) .* cos(y) - 0.52e2 .* x .^ 2 .* exp(-x .^ 2) .* cos(y) + 0.16e2 .* x .^ 4 .* exp(-x .^ 2) .* cos(y) - 0.2e1 .* k .^ 2 .* exp(-x .^ 2) .* cos(y) + 0.4e1 .* k .^ 2 .* x .^ 2 .* exp(-x .^ 2) .* cos(y);
            elseif(dx==0 && dy==2)
                cg = 0.3e1 .* exp(-x .^ 2) .* cos(y) - 0.4e1 .* x .^ 2 .* exp(-x .^ 2) .* cos(y) - k .^ 2 .* exp(-x .^ 2) .* cos(y);
            elseif(dx==1 && dy==1)
                cg = -0.14e2 .* x .* exp(-x .^ 2) .* sin(y) + 0.8e1 .* x .^ 3 .* exp(-x .^ 2) .* sin(y) + 0.2e1 .* k .^ 2 .* x .* exp(-x .^ 2) .* sin(y);
            end
            
        case 7

            
        case 8

            
        case 9

        %%% Plane-wave with k=10, where background has a different k.    
        case 10
            k2 = 10;
            if(dx==0 && dy==0)
                cg = (k^2 - k2^2)*u_exact(x,y,0,0,k,10,false,POL);
            elseif(dx==1 && dy==0)
                cg = (k^2 - k2^2)*u_exact(x,y,1,0,k,10,false,POL);
            elseif(dx==0 && dy==1)
                cg = (k^2 - k2^2)*u_exact(x,y,0,1,k,10,false,POL);
            elseif(dx==2 && dy==0)
                cg = (k^2 - k2^2)*u_exact(x,y,2,0,k,10,false,POL);
            elseif(dx==0 && dy==2)
                cg = (k^2 - k2^2)*u_exact(x,y,0,2,k,10,false,POL);
            elseif(dx==1 && dy==1)
                cg = (k^2 - k2^2)*u_exact(x,y,1,1,k,10,false,POL);
            end
            
        %%% Bump function source.
        case 11
            cg = 0;
            a = 1/2; b = 10;
            if(x.^2+y.^2 >= a^2)
                cg = 0;
            elseif(dx==0 && dy==0)
                cg = exp(-1./(a.^2 - (x.^2+y.^2)));
            elseif((dx==1 && dy==0) || (dx==0 && dy==1))
                cg = -0.2e1 ./ (a ^ 2 - x .^ 2 - y .^ 2) .^ 2 .* x .* exp(-0.1e1 ./ (a .^ 2 - x .^ 2 - y .^ 2));                
            elseif((dx==2 && dy==0) || (dx==0 && dy==2))
                cg = -0.8e1 ./ (a .^ 2 - x .^ 2 - y .^ 2) .^ 3 .* x .^ 2 .* exp(-0.1e1 ./ (a .^ 2 - x .^ 2 - y .^ 2)) - 0.2e1 ./ (a .^ 2 - x .^ 2 - y .^ 2) .^ 2 .* exp(-0.1e1 ./ (a .^ 2 - x .^ 2 - y .^ 2)) + 0.4e1 ./ (a .^ 2 - x .^ 2 - y .^ 2) .^ 4 .* x .^ 2 .* exp(-0.1e1 ./ (a .^ 2 - x .^ 2 - y .^ 2));        
            elseif(dx==1 && dy==1)
                cg = -0.8e1 ./ (a .^ 2 - x .^ 2 - y .^ 2) .^ 3 .* x .* exp(-0.1e1 ./ (a .^ 2 - x .^ 2 - y .^ 2)) .* y + 0.4e1 ./ (a .^ 2 - x .^ 2 - y .^ 2) .^ 4 .* x .* y .* exp(-0.1e1 ./ (a .^ 2 - x .^ 2 - y .^ 2));
            end
            cg = b*cg;
            
        %%% Incident plane-wave
        case 12
%             cg = u_exact(x,y,dx,dy,3,1,false,POL);
            cg = 0;
        case 13
            cg = radiative_source(x,y,k,dx,dy);
        case 14
            cg = radiative_source(x,y,k,dx,dy);
        
        case 20
            cg = 0;
        case 21
            cg = 0;
        case 22
            cg = 0;
        case 23
            cg = 0;
        case 24
            cg = 0;
        case 25
            cg = 0;
        case 26
            cg = 0;
        case 27
            cg = 0;

        case 31
            cg = 0;
        case 32
            cg = 0;
        case 33
            cg = 0;
        case 34
            cg = 0;
        case 35
            cg = 0;
        case 36
            cg = 0;
        case 37
            cg = 0;
        case 38
            cg = 0;
        case 39
            cg = 0;
        case 40
            cg = 0;
        case 50
            cg = 0;
        case 51
            cg = 0;
        case 52
            cg = 0;
        case 53
            A = 1; B = 0.25;
            if(x.^2+y.^2 >= B)
                cg = 0;
            elseif(dx==0 && dy==0)
                cg = A .* exp(-0.1e1 ./ (-x .^ 2 - y .^ 2 + B));
            elseif(dx==1 && dy==0)
                cg = -0.2e1 .* A ./ (-x .^ 2 - y .^ 2 + B) .^ 2 .* x .* exp(-0.1e1 ./ (-x .^ 2 - y .^ 2 + B));
            elseif(dx==0 && dy==1)
                cg = -0.2e1 .* A ./ (-x .^ 2 - y .^ 2 + B) .^ 2 .* y .* exp(-0.1e1 ./ (-x .^ 2 - y .^ 2 + B));
            elseif(dx==1 && dy==1)
                cg = -0.8e1 .* A ./ (-x .^ 2 - y .^ 2 + B) .^ 3 .* x .* exp(-0.1e1 ./ (-x .^ 2 - y .^ 2 + B)) .* y + 0.4e1 .* A ./ (-x .^ 2 - y .^ 2 + B) .^ 4 .* x .* y .* exp(-0.1e1 ./ (-x .^ 2 - y .^ 2 + B));
            elseif(dx==2 && dy==0)
                cg = -0.8e1 .* A ./ (-x .^ 2 - y .^ 2 + B) .^ 3 .* x .^ 2 .* exp(-0.1e1 ./ (-x .^ 2 - y .^ 2 + B)) - 0.2e1 .* A ./ (-x .^ 2 - y .^ 2 + B) .^ 2 .* exp(-0.1e1 ./ (-x .^ 2 - y .^ 2 + B)) + 0.4e1 .* A ./ (-x .^ 2 - y .^ 2 + B) .^ 4 .* x .^ 2 .* exp(-0.1e1 ./ (-x .^ 2 - y .^ 2 + B));
            elseif(dx==0 && dy==2)
                cg = -0.8e1 .* A ./ (-x .^ 2 - y .^ 2 + B) .^ 3 .* y .^ 2 .* exp(-0.1e1 ./ (-x .^ 2 - y .^ 2 + B)) - 0.2e1 .* A ./ (-x .^ 2 - y .^ 2 + B) .^ 2 .* exp(-0.1e1 ./ (-x .^ 2 - y .^ 2 + B)) + 0.4e1 .* A ./ (-x .^ 2 - y .^ 2 + B) .^ 4 .* y .^ 2 .* exp(-0.1e1 ./ (-x .^ 2 - y .^ 2 + B));
            end
        case 54
            cg = 0;
        case 55
            cg = 0;
        case 56
            cg = 0;
    end
end




%       Template for new function
%
%
%
%             if(dx==0 && dy==0)
%                 
%             elseif(dx==1 && dy==0)
%                 
%             elseif(dx==0 && dy==1)
%                 
%             elseif(dx==2 && dy==0)
%                 
%             elseif(dx==0 && dy==2)
%                 
%             end


