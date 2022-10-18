%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Called during init step. This will begin the basic assembly of the
%   important parts of the subdomains. Takes in DOMAIN which is an array of
%   which "blocks" are being used in which locations.
%
%   s is short for start, its called several times so this is a convenience
%   name.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SUBDOMAIN, adj] = subdomain_construction(SUBDOMAIN, DOMAIN, BLOCKS, M, s)
    N = length(SUBDOMAIN);
    e = zeros(N,4);
    [N1, N2] = size(DOMAIN);
    subdomain_counter = 0;
    
    for i=1:N1
        for j=1:N2
            if(DOMAIN(i,j) > 0)
                subdomain_counter = subdomain_counter+1;
                e(subdomain_counter,:) = [s(1)+(j-1)*2, s(2)+(j-1)*2, s(3)-(i-1)*2, s(4)-(i-1)*2];
                SUBDOMAIN(subdomain_counter).edges = e(subdomain_counter,:);
                SUBDOMAIN(subdomain_counter).basis = BLOCKS{DOMAIN(i,j),2};
                SUBDOMAIN(subdomain_counter).type = BLOCKS{DOMAIN(i,j),1};
                SUBDOMAIN(subdomain_counter).params.M = assign_M(SUBDOMAIN(subdomain_counter), M);
                SUBDOMAIN(subdomain_counter).params.k = BLOCKS{DOMAIN(i,j),3};
            end
        end
    end
    adj = subdomain_adjacency(e);
end


%%%---------- Helper Functions ----------%%%
function Mv = assign_M(SUBDOMAIN, M)
    Mv = zeros(1,5);
    for i=1:4
        if(SUBDOMAIN.basis(i) == 'c')
            Mv(i) = M(1);
        elseif(SUBDOMAIN.basis(i) == 'f')
            Mv(i) = M(2);
        end
    end
    if(SUBDOMAIN.type == 1)
        Mv(5) = 0;
    elseif(SUBDOMAIN.type == 2)
        Mv(5) = M(3);
    end    
end
function adj = subdomain_adjacency(S)
    [Ns,~] = size(S);
%     adj = zeros(Ns, 5);
    adj = zeros(Ns,4);
    for i=1:(Ns-1)
        for j=(i+1):Ns
            if(S(i,1) == S(j,2))
                if(S(i,3) == S(j,3))
                    adj(i,1) = j;
                    adj(j,2) = i;
                end
            end
            if(S(i,2) == S(j,1))
                if(S(i,3) == S(j,3))
                    adj(i,2) = j;
                    adj(j,1) = i;
                end
            end
            if(S(i,3) == S(j,4))
                if(S(i,1) == S(j,1))
                    adj(i,3) = j;
                    adj(j,4) = i;
                end
            end
            if(S(i,4) == S(j,3))
                if(S(i,1) == S(j,1))
                    adj(i,4) = j;
                    adj(j,3) = i;
                end
            end
        end
    end


end