%-------------------------------------------------------------------------
%
%   Function: u_exact
%
%   Inputs:     x - x value for evaluation
%               y - y value for evaluation
%              dx - number of x derivatives
%              dy - number of y derivatives
%               r - which test function to use
%                      0 - custom
%                      1 ... 8 - preset, see 'test_function.m'
%
%   Output:     u - indicated function evaluated at (x,y)
%
%   Description:    This function provides a central location from which
%                   from which to obtain the exact solution, no matter
%                   which case is being tested.
%   
%-------------------------------------------------------------------------

function u = u_exact(x,y,dx,dy,k,r,REAL,POL)
    
    if(exist('POL', 'var') == 0)
        POL = false;
    end
    
    if(~POL)
        u = test_function(x,y,dx,dy,k,r,false);
    else
        rad = x; t = y;
        [x,y] = pol2cart(t,rad);
%         u = 0;
        a = 0; b = 0;
        if(dx==0 && dy==0)
            u = test_function(x+a,y+b,0,0,k,r,false);
        elseif(dx==1 && dy==0)
            % first "r" derivative
            u = test_function(x+a,y+b,1,0,k,r,false) .* cos(t) + ...
                 test_function(x+a,y+b,0,1,k,r,false) .* sin(t);
        elseif(dx==2 && dy==0)
            % second "r" derivative
            u = test_function(x+a,y+b,2,0,k,r,false) .* (cos(t)).^2 + ...
                 test_function(x+a,y+b,0,2,k,r,false) .* (sin(t)).^2 + ...
                 2*test_function(x+a,y+b,1,1,k,r,false) .* cos(t) .* sin(t);
        elseif(dx==0 && dy==1)
            u = test_function(x+a,y+b,1,0,k,r,false) .* (-rad.*sin(t)) + ...
                 test_function(x+a,y+b,0,1,k,r,false) .* (rad.*cos(t));
        elseif(dx==0 && dy==2)
            u = test_function(x+a,y+b,2,0,k,r,false) .* (-rad.*sin(t)).^2 + ...
                 test_function(x+a,y+b,0,2,k,r,false) .* (rad.*cos(t)).^2 + ...
                 test_function(x+a,y+b,1,1,k,r,false) .* 2.*(-rad.^2.*cos(t).*sin(t)) + ...
                 test_function(x+a,y+b,1,0,k,r,false) .* (-rad.*cos(t)) + ...
                 test_function(x+a,y+b,0,1,k,r,false) .* (-rad.*sin(t));
        end
    end
    %%% This is pretty much a moot point now.
    if(REAL)
        u = real(u);
    end
end

%-------------------------------------------------------------------------
%
%   Function: test_function
%
%   Inputs:     x - x value for evaluation
%               y - y value for evaluation
%              dx - number of x derivatives
%              dy - number of y derivatives
%               r - which test function to use
%                      1 ... 8 - preset, see 'test_function.m'
%
%   Output:     u - indicated function evaluated at (x,y)
%
%   Description:    Houses all pre-written test cases. The function itself,
%                   as well as the first derivative in each direction is 
%                   provided, so that neumann boundary data is viable.
%   
%-------------------------------------------------------------------------


function u = test_function(x,y,dx,dy,k,r,POL)
    switch(r)
        
        %%% Simplest test cases (plane-waves and simple inhomogeneous)
        case 1
            k = k(1);
            theta = 0;
            c = cos(theta);
            s = sin(theta);
            if(dx==0 && dy==0)
                u = exp(1i .* k .* (c*x + s*y));
            elseif(dx==1 && dy==0)
                u = 1i .* k .* c .* exp(1i .* k .* (c*x + s*y));
            elseif(dx==0 && dy==1)
                u = 1i .* k .* s .* exp(1i .* k .* (c*x + s*y));
            elseif(dx==2 && dy==0)
                u = -k^2 .* c^2 .* exp(1i .* k .* (c*x + s*y));
            elseif(dx==0 && dy==2)
                u = -k^2 .* s^2 .* exp(1i .* k .* (c*x + s*y));
            elseif(dx==1 && dy==1)
                u = -k^2 .* s .* c .* exp(1i .* k .* (c*x + s*y));
            end
        case 2
            k = k(1);
            if(POL)
                r = x; t = y;
                if(dx==0 && dy==0)
                    u = exp(1i.*k.*r.*sin(t));
                elseif(dx==1 && dy==0)
                    u = sin(t) .* 1i .* k .* exp(1i.*k.*r.*sin(t));
                elseif(dx==0 && dy==1)
                    u = 1i*k*r.*cos(t) .* exp(1i.*k.*r.*sin(t));
                end
            else
                if(dx==0 && dy==0)
                    u = exp(1i.*k.*y);
                elseif(dx==1 && dy==0)
                    u = 0;
                elseif(dx==0 && dy==1)
                    u = 1i.*k.*exp(1i.*k.*y);
                end
            end
            
        case 3
            k = k(1);
            theta = pi/4;
            c = cos(theta); s = sin(theta);
            if(dx==0 && dy==0)
                u = exp(1i .* k .* (c*x + s*y));
            elseif(dx==1 && dy==0)
                u = 1i .* k .* c .* exp(1i .* k .* (c*x + s*y));
            elseif(dx==0 && dy==1)
                u = 1i .* k .* s .* exp(1i .* k .* (c*x + s*y));
            elseif(dx==2 && dy==0)
                u = -k^2 .* c^2 .* exp(1i .* k .* (c*x + s*y));
            elseif(dx==0 && dy==2)
                u = -k^2 .* s^2 .* exp(1i .* k .* (c*x + s*y));
            elseif(dx==1 && dy==1)
                u = -k^2 .* s .* c .* exp(1i .* k .* (c*x + s*y));
            end            
            
        case 4
            if(POL)
                r = x; t = y;
                if(r >= 1)
                    u = 0;
                elseif(dx==0 && dy==0)
                    u = exp(-0.1e1 ./ (1 - (r.^2)));
                elseif(dx==1 && dy==0)
                    u = -0.2e1 / ((-(r .^ 2) + 1) .^ 2) .* r .* exp(-0.1e1 / (-(r .^ 2) + 1));
                elseif(dx==0 && dy==1)
                    u = 0;
                end
            else
                if(x^2+y^2>=1)
                    u = 0;
                elseif(dx==0 && dy==0)
                    u = exp(-0.1e1 ./ (1 - y .^ 2 - x .^ 2));
                elseif(dx==1 && dy==0)
                    u = -0.2e1 ./ ((-x .^ 2 - y .^ 2 + 1) .^ 2) .* x .* exp(-0.1e1 ./ (-x .^ 2 - y .^ 2 + 1));
                elseif(dx==0 && dy==1)
                    u = -0.2e1 ./ ((-x .^ 2 - y .^ 2 + 1) .^ 2) .* y .* exp(-0.1e1 ./ (-x .^ 2 - y .^ 2 + 1));
                end
            end
            
        case 5
            wave1 = 1;
            wave2 = 3;
            
            u = u_exact(x,y,dx,dy,k,wave1,false,false) + u_exact(x,y,dx,dy,k,wave2,false,false);

        case 6
            if(POL)
                r = x; t = y;
                if(dx==0 && dy==0)
                    u = exp(-r .^ 2 .* cos(t) .^ 2) .* cos(r .* sin(t));
                elseif(dx==1 && dy==0)
                    u = -0.2e1 * r .* cos(t) .^ 2 .* exp(-r .^ 2 .* cos(t) .^ 2) .* cos(r .* sin(t)) - exp(-r .^ 2 .* cos(t) .^ 2) .* sin(t) .* sin(r .* sin(t));
                end
            else
                if(dx==0 && dy==0)
                    u = exp(-x.^2) .* cos(y);
                elseif(dx==1 && dy==0)
                    u = -0.2e1 .* x .* exp(-x .^ 2) .* cos(y);
                elseif(dx==0 && dy==1)
                    u = -exp(-x .^ 2) .* sin(y);
                end
            end
            
        case 7

            
        case 8

            
        case 9

            
        %%% Two-wavenumber case (simple waves)
        case 10
            k=10;
            u = u_exact(x,y,dx,dy,k,5,false,POL);
        %%% Self-convergence with a bump source, set boundary condition to zero
        case 11
            u = 0;
        %%% Incident wave, this zero BC is for an ABC
        case 12
            u = u_exact(x,y,dx,dy,k,1,false,POL);
        %%% Used to test exact radiation solution
        case 13
            u = radiative_solution(x,y,k,dx,dy);
        %%% Used to approximate radiation solution
        case 14
            u = 0; % RHS of ABC
            
        %%% Bessel function cases
        case 20
            if(dx==0 && dy==0)
                u = sin(k(1)*x);
            elseif(dx==1 && dy==0)
                u = k(1) * cos(k(1)*x);
            else
                u = 0;
            end
        case 21
            if(dx==0 && dy==0)
                u = sin(k(1)*y);
            elseif(dx==0 && dy==1)
                u = k(1) * cos(k(1)*y);
            else
                u = 0;
            end
        
            


        case 31
            a = [1, zeros(1,9)];
            u = bessel_eval(x,y,2./3,k,a);
        case 32
            a = [0, 1, zeros(1,8)];
            u = bessel_eval(x,y,2./3,k,a);
        case 33
            a = [0, 0, 1, zeros(1,7)];
            u = bessel_eval(x,y,2./3,k,a);
        case 34
            a = [0 0 0 1 0 0 0 0 0 0];
            u = bessel_eval(x,y,2./3,k,a);
        case 35
            a = [0 0 0 0 1 0 0 0 0 0];
            u = bessel_eval(x,y,2./3,k,a);
        case 36
            a = [zeros(1,5), 1, zeros(1,4)];
            u = bessel_eval(x,y,2./3,k,a);
        case 37
            a = [zeros(1,6), 1, zeros(1,3)];
            u = bessel_eval(x,y,2./3,k,a);
        case 38
            a = [zeros(1,7), 1, 0 0];
            u = bessel_eval(x,y,2./3,k,a);
        case 39
            a = [zeros(1,8), 1, 0];
            u = bessel_eval(x,y,2./3,k,a);
        case 40
            a = [zeros(1,9), 1];
            u = bessel_eval(x,y,2./3,k,a);
            
        case 50 % case for two wavenumbers, straight horizontal
            R = (k(1) - k(2)) ./ (k(1) + k(2));
            if(x<=0)
                u = exp(1i.*k(1).*(-x)) + ...
                    R .* exp(1i.*k(1).*x);
            else
                T = 1 + R;
                u = T .* exp(1i.*k(2).*(-x));
            end
            
        case 51
            t = pi./4;
            tt = asin((k(1)./k(2)) .* sin(t));
            Kr = (cos(tt).*k(2)) ./ (cos(t).*k(1));
            if(dx==0 && dy==0)
                if(x<=0)
                    R = (1 - Kr) ./ (1 + Kr);
                    u = exp(1i.*k(1).*(-cos(t).*x - sin(t).*y)) + ...
                        R .* exp(1i.*k(1).*(cos(t).*x - sin(t).*y));
                else
                    T = 2 ./ (1+Kr);
                    u = T .* exp(1i.*k(2).*(-cos(tt).*x - sin(tt).*y));
                end
            elseif(dx==1 && dy==0)
                if(x<=0)
                    R = (1 - Kr) ./ (1 + Kr);
                    u = exp(1i.*k(1).*(-cos(t).*x - sin(t).*y)) .* (-1i.*k(1).*cos(t)) + ...
                        R .* exp(1i.*k(1).*(cos(t).*x - sin(t).*y)) .* (1i.*k(1).*cos(t));
                else
                    T = 2 ./ (1+Kr);
                    u = T .* exp(1i.*k(2).*(-cos(tt).*x - sin(tt).*y)) .* (-1i.*k(2).*cos(tt));
                end
            elseif(dx==0 && dy==1)
                if(x<=0)
                    R = (1 - Kr) ./ (1 + Kr);
                    u = exp(1i.*k(1).*(-cos(t).*x - sin(t).*y)) .* (-1i.*k(1).*sin(t)) + ...
                        R .* exp(1i.*k(1).*(cos(t).*x - sin(t).*y)) .* (1i.*k(1).*sin(t));
                else
                    T = 2 ./ (1+Kr);
                    u = T .* exp(1i.*k(2).*(-cos(tt).*x - sin(tt).*y)) .* (-1i.*k(1).*sin(tt));
                end
            end
            
        case 52
            R = 1;
            MM = 1;
            [th, r] = cart2pol(x,y);
        %    u = -(besselj(0,k(1).*R) ./ besselh(0,k(1).*R)) .* besselh(0,k(1).*r);
            u = 0;
            for m = 1:MM
                u = u - 2.*(1i).^m .* (besselj(m,k(1).*R) ./ besselh(m,k(1).*R)) .* besselh(m,k(1).*r) .* cos(m.*th);
            end
            
        case 53 
            u = 0;
            
        case 54 % k = [3, 5, 13, 20], x = [-1, 1, 3, 5, 7]
              c= [1.000000000000000 + 0.000000000000000i;
                 -0.824052738201097 - 0.566513093108425i;
                 -0.421034409081433 - 0.547895732182119i;
                  0.657344663619833 - 0.212973472959569i;
                  0.088736651005524 - 0.257493668189648i;
                  0.072749854182055 - 0.262458836974347i;
                 -0.021871117507080 + 0.213663370454081i;
                 -0.103020142610685 + 0.188460159904634i];
            if(x<=1)
                u = c(1).*exp(1i.*3.*x) + c(2).*exp(-1i.*3.*x);
            elseif(x<=3)
                u = c(3).*exp(1i.*5.*x) + c(4).*exp(-1i.*5.*x);
            elseif(x<=5)
                u = c(5).*exp(1i.*13.*x) + c(6).*exp(-1i.*13.*x);
            else
                u = c(7).*exp(1i.*20.*x) + c(8).*exp(-1i.*20.*x);
            end
            
        case 55
            c = two_way_coefficients(k, [-1, 1, 3]);
            if(x<=1)
                u = c(1).*exp(1i.*k(1).*x) + c(2).*exp(-1i.*k(1).*x);
            else
                u = c(3).*exp(1i.*k(2).*x);
            end
        
        case 56
            if(x == -1)
                u = exp(1i .* (k(1)./sqrt(2)) .* (x + y));
            else
                u = 0;
            end
    end
end

%    Template for a new function
%
%             if(dx==0 && dy==0)
%                 
%             elseif(dx==1 && dy==0)
%                 
%             elseif(dx==0 && dy==1)
%                 
%             end
