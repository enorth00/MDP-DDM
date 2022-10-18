% Reflection errors and their relative percents at the exterior wall
clear all;

%% Set everything for the first one. 1x1 case

OPTIONS.T = 1;

Mc = 30;

k1 = 10;
k2 = 0;

Mf=31; Mr=35; M = [Mc, Mf, Mr];
TYPES = cell(1,3);
TYPES(1,:) = [{1}, {'cccc'}, {[k1, 0]}];
DOMAIN = 1*ones(1); 
start = [-1, 1, -1, 1];

OPTIONS.Ns = numel(DOMAIN);
OPTIONS.TYPES = TYPES;
BLOCK = allocate_blocks([1, 2]);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);


OPTIONS.GRIDS = [6:9];
OPTIONS.FUNC = 14;
OPTIONS.PARALLEL = true;
OPTIONS.REAL = false;
OPTIONS.RADIATION = 2;
OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);
OPTIONS.FOURIER_AUX = [-3*pi/8, 3*pi/8];

run('DRIVER_exact');

errors1 = errors;

u1 = @(y) u_exact(-1,y, 0, 0, k1, 13, false, false);
u2 = @(y) u_exact(1,y, 0, 0, k1, 13, false, false);
u3 = @(x) u_exact(x,-1, 0, 0, k1, 13, false, false);
u4 = @(x) u_exact(x,1, 0, 0, k1, 13, false, false);
cheb1 = real(chebfun(u1, [-1, 1]));
cheb2 = real(chebfun(u2, [-1, 1]));
cheb3 = real(chebfun(u3, [-1, 1]));
cheb4 = real(chebfun(u4, [-1, 1]));

bound_max = max([max(cheb1), max(cheb2), max(cheb3), max(cheb4)]);

percents1 = errors1 / bound_max;





%% 3x3 case

DOMAIN = 1*ones(3);
start = [-3, -1, 1, 3];
OPTIONS.Ns = numel(DOMAIN);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);
OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);

run('DRIVER_exact');

errors2 = errors;

u1 = @(y) u_exact(-3,y, 0, 0, k1, 13, false, false);
u2 = @(y) u_exact(3,y, 0, 0, k1, 13, false, false);
u3 = @(x) u_exact(x,-3, 0, 0, k1, 13, false, false);
u4 = @(x) u_exact(x,3, 0, 0, k1, 13, false, false);
cheb1 = real(chebfun(u1, [-3, 3]));
cheb2 = real(chebfun(u2, [-3, 3]));
cheb3 = real(chebfun(u3, [-3, 3]));
cheb4 = real(chebfun(u4, [-3, 3]));

bound_max = max([max(cheb1), max(cheb2), max(cheb3), max(cheb4)]);

percents2 = errors2 / bound_max;

%% 5x5 case

DOMAIN = 1*ones(5);
start = [-5, -3, 3, 5];
OPTIONS.Ns = numel(DOMAIN);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);
OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);

run('DRIVER_exact');

errors3 = errors;

u1 = @(y) u_exact(-5,y, 0, 0, k1, 13, false, false);
u2 = @(y) u_exact(5,y, 0, 0, k1, 13, false, false);
u3 = @(x) u_exact(x,-5, 0, 0, k1, 13, false, false);
u4 = @(x) u_exact(x,5, 0, 0, k1, 13, false, false);
cheb1 = real(chebfun(u1, [-5, 5]));
cheb2 = real(chebfun(u2, [-5, 5]));
cheb3 = real(chebfun(u3, [-5, 5]));
cheb4 = real(chebfun(u4, [-5, 5]));

bound_max = max([max(cheb1), max(cheb2), max(cheb3), max(cheb4)]);

percents3 = errors3 / bound_max;


output = [errors1, percents1, errors2, percents2, errors3, percents3];

disp(output);


%% output a latex table
n = (2 .^ [6:9])';
output = [n, errors1, 100*percents1, errors2, 100*percents2, errors3, 100*percents3];

formatspec = "%i & %.2e & %.2f\\%% & %.2e & %.2f\\%% & %.2e & %.2f\\%% \\\\\n";

fprintf(formatspec, output');
