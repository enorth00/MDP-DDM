function d = basiscoeffs(f,M,basis)

    if(basis == 'c')
        d = chebcoeffs(f, M);
    elseif(basis == 'f')
%         d = trigcoeffs(f, M);
        d = f.coeffs;
        F = length(f);
        if(M>F)
            if(mod(M-F,2)==0)
                d = [zeros((M-F)/2,1); d; zeros((M-F)/2,1)];
            else
                disp('Length of f was weird in basiscoeffs');
            end
        elseif(M<F)
            if(mod(F-M,2)==0)
                d = d((F+1)/2 - (M-1)/2 : (F+1)/2 + (M-1)/2);
            else
                disp('Length of f was weird in basiscoeffs');
            end
        end
    end

end
