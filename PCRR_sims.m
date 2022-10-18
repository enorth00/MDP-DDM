% PCRR_sims.m

clearvars;
% close all;

domain = 'PCRR_sims_init';

run(['inits/', domain]);
% [OPTIONS, BLOCK] = set_GlobalOptions(OPTIONS);

% Set the discretizations, test functions, parallel, and BCs.
OPTIONS.GRIDS = [7];
OPTIONS.FUNC = 14;
OPTIONS.PARALLEL = true;
OPTIONS.REAL = false;
OPTIONS.RADIATION = 2;
OPTIONS.FOURIER_AUX = [-3*pi/8, 3*pi/8];
OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);

run('DRIVER_self');

save('last_run.mat', 'SOLUTION', 'BLOCK', 'OPTIONS', 'SUBDOMAIN');
