% MAKE_new_Q.m
clear;

domain = 'empty_base';
% domain = 'rod_base';

run(['inits/', domain]);

% Kvals = [10]'; K = [Kvals, zeros(length(Kvals),1)];
% K = [1.1, 1.1*3.48];
% K = [1.1, 0];
K = [3, 0];
GRIDS = [5:9];
% M = [100, 51];
M = [30, 0];

T = '';
switch(SUBDOMAIN.type)
    case 1
        T = 'empty';
    case 2
        T = 'rod';
end

OPTIONS.PARALLEL = false;
OPTIONS.REAL = false;

% OPTIONS.FOURIER_AUX = [-pi, pi];
% OPTIONS.FOURIER_AUX = [-7*pi/8, 7*pi/8];
% OPTIONS.FOURIER_AUX = [-pi/2, pi/2];
OPTIONS.FOURIER_AUX = [-3*pi/8, 3*pi/8];
% OPTIONS.FOURIER_AUX = [-1.1, 1.1];
% OPTIONS.FOURIER_AUX = [-1, 1];

disp('Starting Construction...');

[K_length,~] = size(K);
SUBDOMAIN.params.M = [M(1), M(1), M(1), M(1), M(2)];

for j = 1:K_length
SUBDOMAIN.params.k = K(j,:);
fprintf("Wavenumber = %i\n", K(j,1));

for i = 1:length(GRIDS)
    N = 2^(GRIDS(i))+1;
    fprintf('\nGrid %i ... ', N-1);
    
    
%     k = K;
%     SUBDOMAIN.params.k = K;
    BLOCK(SUBDOMAIN.type).N = [N,N];
    BLOCK(SUBDOMAIN.type).grids = SET_grids(BLOCK(SUBDOMAIN.type));
    
    bf = SET_basis_functions(M, 'chebyshev', OPTIONS);
    [Q, Qrod] = construct_Q(BLOCK(SUBDOMAIN.type), SUBDOMAIN, bf, OPTIONS);
    Q_store.Q = Q; Q_store.Qrod = Qrod;
    save(['Q_mats/' T '/chebyshev/Q_' num2str(K(j,1)) '_' num2str(K(j,2)) '_' num2str(N-1) '.mat'], 'Q_store');
    
%     bf = SET_basis_functions(M, 'fourier', OPTIONS);
%     [Q, Qrod] = construct_Q(BLOCK(SUBDOMAIN.type), SUBDOMAIN, bf, OPTIONS);
%     save(['Q_mats/' T '/fourier/Q_' num2str(K(j,1)) '_' num2str(K(j,2)) '_' num2str(N-1) '.mat'], 'Q');

    fprintf('(%i, %i)', K(j,:));
    
    disp('');
    disp('Done');
end
end
