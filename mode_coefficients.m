function [a,b] = mode_coefficients(N, wavenumber, OPTIONS)

    CONDITION = OPTIONS.RADIATION;
    D = OPTIONS.FOURIER_AUX;
    L = D(2) - D(1);

    % 1 = Sommerfeld
    % 2 = BTG 2nd-order
    % 3 = Fourier mode-canceling
    
    a = zeros(2*N+1,1);
    b = ones(2*N+1,1);
    
    if(CONDITION == 1) %Sommerfeld
        a = 1i*wavenumber*ones(2*N+1,1);
    else
        for mode=0:N % mode counter for indexing
            
            mode2 = (2*pi/L) * mode;
            if(CONDITION == 5) % Fourier BTG 2nd-order
                a(N+1-mode) = (1/2)*mode2^2 - wavenumber^2;
                a(N+1+mode) = a(N+1-mode);
                b(N+1-mode) = -1i*wavenumber;
                b(N+1+mode) = -1i*wavenumber;
                
            elseif(CONDITION == 6) %Fourier mode-canceling
%                 if(mode2 > wavenumber)
%                     mode2=0;
%                 end
                if(abs(mode2^2 - wavenumber^2) < 1e-8)
                    alpha = 1;
                    beta(N+1-mode) = 0;
                    beta(N+1+mode) = 0;
                elseif(mode2 >= wavenumber) % evanescent case
%                     alpha = -1i*wavenumber;
                    alpha = sqrt(mode2^2 - wavenumber^2);
                else % oscillating case
                    alpha = -1i*sqrt(wavenumber^2 - mode2^2);
                end
                a(N+1-mode) = alpha;
                a(N+1+mode) = alpha;                
            else
                disp('CONDITION is invalid');
            end
        end
    end

end
