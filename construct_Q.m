%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Function: construct_Q
%
%   Input:  u - (M+1)x(N+1) matrix 
%      params - struct containing at least: N, h, k
%
%   Output: g - (M+1)x(N+1) matrix of stenciled values.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q, Qrod] = construct_Q(BLOCK, SUBDOMAIN, bf, OPTIONS)
    R = OPTIONS.REAL;  
    Ng = BLOCK.grids.Ngamma;
    M = SUBDOMAIN.params.M(1:4);
    if(BLOCK.radius>0)
        Mr = SUBDOMAIN.params.M(5);
    else
        Mr = 0;
    end
    
    NN = length(Ng(Ng));
    
    % total basis functions
    tbf0 = 2*sum(SUBDOMAIN.params.M);
    
    extensions0 = extend_square_basisfuncs(bf.b0, 0, bf.b2, 0, bf.b4,...
                                           bf.r, BLOCK, SUBDOMAIN, ...
                                           SUBDOMAIN.params.k(1), 0);
    extensions1 = extend_square_basisfuncs(0, bf.b0, 0, bf.b2, 0, ...
                                           bf.r, BLOCK, SUBDOMAIN, ...
                                           SUBDOMAIN.params.k(1), 1);
    extensions = [extensions0, extensions1];
    if(Mr>0)
        extrod0full = extend_square_basisfuncs(bf.b0, 0, bf.b2, 0, bf.b4, ...
                                               bf.r, BLOCK, SUBDOMAIN, ...
                                               SUBDOMAIN.params.k(2), 0);
        extrod1full = extend_square_basisfuncs(0, bf.b0, 0, bf.b2, 0, ...
                                               bf.r, BLOCK, SUBDOMAIN, ...
                                               SUBDOMAIN.params.k(2), 1);
        
        extrod0 = extrod0full(BLOCK.grids.coords(:,3)==5,sum(M)+1:sum(M)+Mr);
        extrod1 = extrod1full(BLOCK.grids.coords(:,3)==5,sum(M)+1:sum(M)+Mr);

        extrod = [extrod0, extrod1];
        proj_rod = zeros(size(extrod));
    end     
    fprintf('Extensions complete!\n');

%     extensions_short = [extensions0(:,1:M), extensions0(:,2*M+1:3*M), extensions1(:,1:M), extensions1(:,2*M+1:3*M)];
%     [~, tbf0] = size(extensions_short);
    projections = zeros(NN,tbf0);
    if(true)
    if(OPTIONS.PARALLEL)
        parfor K = 1:tbf0
            w = zeros(size(Ng));
    %         w(Ng) = extensions_short(:,K);
            w(Ng) = extensions(:,K);
            Pw = N_plus_potential_Q(w, BLOCK, SUBDOMAIN, R);
            projections(:,K) = Pw(Ng);
        end
        
        if(Mr>0)
            [Nrod, tbfrod] = size(extrod);
            rNg = BLOCK.grids.rodNgamma;
            for K = 1:tbfrod
                w = zeros(size(rNg));
        %         w(rNg) = extensions_short(:,K);
                w(rNg) = extrod(:,K);
                Pw = N_plus_potential_rod(w, BLOCK, SUBDOMAIN, R);
                proj_rod(:,K) = Pw(rNg);
            end
        end
    else
        textprogressbar('Square projections: ');
        for K = 1:tbf0
            w = zeros(size(Ng));
    %         w(Ng) = extensions_short(:,K);
            w(Ng) = extensions(:,K);
            Pw = N_plus_potential_Q(w, BLOCK, SUBDOMAIN, R);
            projections(:,K) = Pw(Ng);
            textprogressbar(100*(K/tbf0));
        end
        textprogressbar('done');
        
        if(Mr>0)
            [Nrod, tbfrod] = size(extrod);
            textprogressbar('Rod projections: ');
            rNg = BLOCK.grids.rodNgamma;
            for K = 1:tbfrod
                w = zeros(size(rNg));
        %         w(rNg) = extensions_short(:,K);
                w(rNg) = extrod(:,K);
                Pw = N_plus_potential_rod(w, BLOCK, SUBDOMAIN, R);
                proj_rod(:,K) = Pw(rNg);
                textprogressbar(100*(K/tbfrod));
            end
            textprogressbar('done');
        end
    end
    end
%     s = BLOCK.grids.sides;
%     Qt = projections - extensions_short;
%     Qt2 = rotate_functions(Qt(:,1:M),2,s); Qt4 = rotate_functions(Qt(:,M+1:2*M), 2, s);
%     Qt12 = rotate_functions(Qt(:,2*M+1:3*M), 2, s); Qt14 = rotate_functions(Qt(:,3*M+1:4*M), 2, s);
%     
%     Q = [Qt(:,1:M), Qt2, Qt(:,M+1:2*M), Qt4, Qt(:,2*M+1:3*M), Qt12, Qt(:,3*M+1:4*M), Qt14];


    Q = projections - extensions;
%     Q = extensions;
    
    
    if(Mr>0)
        Qrod = proj_rod - extrod;
    else
        Qrod = [];
    end
end
    
    
    