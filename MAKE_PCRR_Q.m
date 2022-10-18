% MAKE_PCRR_Q.m
clear;

K1 = 1;
GRIDS = [6:9];
M1 = 70;
M2 = 51;

% K2 = K1*3.48;
% K2 = K1;
K2 = 3;

OPTIONS.PARALLEL = false;
OPTIONS.REAL = false;
OPTIONS.FOURIER_AUX = [-3*pi/8, 3*pi/8];

disp('Starting Construction...');

for i=1:length(GRIDS)
    N = 2^(GRIDS(i)) + 1;
    fprintf('\nGrid %i ... ', N-1);

    run('inits/empty_base');
    SUBDOMAIN.params.k = [K1, 0];
    M = [M1*ones(1,4), 0];
    SUBDOMAIN.params.M = M;
    BLOCK(SUBDOMAIN.type).N = [N,N];
    BLOCK(SUBDOMAIN.type).grids = SET_grids(BLOCK(SUBDOMAIN.type));
    bf = SET_basis_functions([M1, 0], 'chebyshev', OPTIONS);
    [Q, Qrod] = construct_Q(BLOCK(SUBDOMAIN.type), SUBDOMAIN, bf, OPTIONS);
    Q_store.Q = Q; Q_store.Qrod = Qrod;
    save(['Q_mats/empty/chebyshev/Q_' num2str(K1) '_0_' num2str(N-1) '.mat'], 'Q_store');
    
    run('inits/rod_base');
    SUBDOMAIN.params.k = [K1, K2];
    M = [M1*ones(1,4), M2];
    SUBDOMAIN.params.M = M;
    BLOCK(SUBDOMAIN.type).N = [N,N];
    BLOCK(SUBDOMAIN.type).grids = SET_grids(BLOCK(SUBDOMAIN.type));
    bf = SET_basis_functions([M1, M2], 'chebyshev', OPTIONS);
    [Q, Qrod] = construct_Q(BLOCK(SUBDOMAIN.type), SUBDOMAIN, bf, OPTIONS);
    Q_store.Q = Q; Q_store.Qrod = Qrod;
    save(['Q_mats/rod/chebyshev/Q_' num2str(K1) '_' num2str(K2) '_' num2str(N-1) '.mat'], 'Q_store');
    
end

disp('Done');


