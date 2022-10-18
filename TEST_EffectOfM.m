% TEST_EffectOfM.m

clearvars;
% close all;

domain = 'M_test';

% Set the discretizations, test functions, parallel, and BCs.
OPTIONS.GRIDS = [7];
OPTIONS.FUNC = 14;
OPTIONS.PARALLEL = false;
OPTIONS.REAL = false;
% OPTIONS.RADIATION = 2;

% OPTIONS.FOURIER_AUX = [-pi, pi];
% OPTIONS.FOURIER_AUX = [-7*pi/8, 7*pi/8];
% OPTIONS.FOURIER_AUX = [-pi/2, pi/2];
OPTIONS.FOURIER_AUX = [-3*pi/8, 3*pi/8];
% OPTIONS.FOURIER_AUX = [-1.1, 1.1];
% OPTIONS.FOURIER_AUX = [-1,1];

MF = [15];
RC = [2];
e_vec = zeros(length(MF),4);
for MCASE = 1:length(MF)
    Mc = MF(MCASE);

    for RAD_CASE = RC
%    	Mf = MF;
    	run('inits/M_test');
    	OPTIONS.RADIATION = RAD_CASE;
    	OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
        OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);
       fprintf("M = %i", Mc);
%        fprintf(" --- Radiation Condition = %i\n", OPTIONS.RADIATION);
        run('DRIVER_self');
        e_vec(MCASE,RAD_CASE) = errors(end);
    end
%     pause;
end

% save('last_run.mat', 'SOLUTION', 'BLOCK', 'OPTIONS', 'SUBDOMAIN');

% run('PLOT_server');
% run('PLOT_error');
