function F = apply_B(f,stencil)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Function: apply_B
%   This function applies the right-hand side, 5-node
%   stencil that maps K_0 to G_0 for use in the solver.
%
%   Inputs:     f - M x N matrix of source values
%   
%   Outputs:    F - M x N stenciled source values (with zeros around
%                   outside)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(nargin > 1)
        option = stencil;
    else
        option = true;
    end
    
    [M,N] = size(f);
    F = zeros(M,N);
    
    if(option)
        
        for m = 2:M-1
            for n = 2:N-1
                F(m,n) = 8*f(m,n) + f(m,n+1) + f(m+1,n) ...
                    + f(m,n-1) + f(m-1,n);
            end
        end
        F = F/12;
    else
        F = f;
    end
        

end