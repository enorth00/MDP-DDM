% TEST_exact.m
% 
% Quick test script to initialize what I need for the scatter rod 
% base subdomain
clearvars;
close all;

% pick a DOMAIN to run

% domain = 'empty_base';
% domain = 'empty_base_2';
% domain = 'empty_square_2x2';
% domain = 'empty_square_NxN';
domain = 'playground';

% domain = 'rod_base';
% domain = 'rod_base_2';
% domain = 'rod_square_2x2';

% domain = 'rod_and_empty';
% domain = 'tunnel_3x3';

run(['inits/', domain]);
% [OPTIONS, BLOCK] = set_GlobalOptions(OPTIONS);

% Set the discretizations, test functions, parallel, and BCs.
OPTIONS.GRIDS = [6];
OPTIONS.FUNC = 13;
OPTIONS.PARALLEL = false;
OPTIONS.REAL = false;
OPTIONS.RADIATION = 0;
OPTIONS.FOURIER_AUX = [-3*pi/8, 3*pi/8];

OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);
% OPTIONS.BOUNDARY.beta = 1i*k1*ones(OPTIONS.Ns,4);

% OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
% OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);
 
run('DRIVER_exact');

% run('PLOT_solution');
% run('PLOT_error');
