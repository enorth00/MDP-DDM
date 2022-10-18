%% Case 1. Empty and rod, centered on empty.
clear all;

% k=5

OPTIONS.T = 2;

Mc = 20;
Mr = 31;

k1 = 5;
k2 = 5;

Mf=31;  M = [Mc, Mf, Mr];
TYPES = cell(1,3);
TYPES(1,:) = [{1}, {'cccc'}, {[k1, 0]}];
TYPES(2,:) = [{2}, {'cccc'}, {[k1, k2]}];
DOMAIN = [1, 2]; 
start = [-1, 1, -1, 1];

OPTIONS.Ns = numel(DOMAIN);
OPTIONS.TYPES = TYPES;
BLOCK = allocate_blocks([1, 2]);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);


OPTIONS.GRIDS = [6];
OPTIONS.FUNC = 14;
OPTIONS.PARALLEL = true;
OPTIONS.REAL = false;
OPTIONS.RADIATION = 2;
OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);
OPTIONS.FOURIER_AUX = [-3*pi/8, 3*pi/8];

run('DRIVER_exact');

errors1 = zeros(5,1);
errors1(1) = errors;


Mc = 20;
M = [Mc, Mf, Mr];
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);
OPTIONS.GRIDS = [7:10];
run('DRIVER_exact');

errors1(2:end) = errors;

ratios1 = log2(errors1(1:end-1) ./ errors1(2:end));

%% Case 2

% k = 10

Mc = 20;
Mr = 31;

k1 = 10;
k2 = 10;

M = [Mc, Mf, Mr];
TYPES = cell(1,3);
TYPES(1,:) = [{1}, {'cccc'}, {[k1, 0]}];
TYPES(2,:) = [{2}, {'cccc'}, {[k1, k2]}];
DOMAIN = [1, 2]; 
start = [-1, 1, -1, 1];

OPTIONS.Ns = numel(DOMAIN);
OPTIONS.TYPES = TYPES;
BLOCK = allocate_blocks([1, 2]);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);

OPTIONS.GRIDS = [6];
OPTIONS.FUNC = 13;
OPTIONS.PARALLEL = true;
OPTIONS.REAL = false;
OPTIONS.RADIATION = 0;
OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);
OPTIONS.FOURIER_AUX = [-3*pi/8, 3*pi/8];

run('DRIVER_exact');

errors2 = zeros(5,1);
errors2(1) = errors;

Mc = 30;
M = [Mc, Mf, Mr];
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);
OPTIONS.GRIDS = [7:10];

run('DRIVER_exact');

errors2(2:end) = errors;

ratios2 = log2(errors2(1:end-1) ./ errors2(2:end));


%% Case 3

% k = 20

Mc = 30;
Mr = 41;

k1 = 20;
k2 = 20;

M = [Mc, Mf, Mr];
TYPES = cell(1,3);
TYPES(1,:) = [{1}, {'cccc'}, {[k1, 0]}];
TYPES(2,:) = [{2}, {'cccc'}, {[k1, k2]}];
DOMAIN = [1, 2]; 
start = [-1, 1, -1, 1];

OPTIONS.Ns = numel(DOMAIN);
OPTIONS.TYPES = TYPES;
BLOCK = allocate_blocks([1, 2]);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);

OPTIONS.GRIDS = [6];
OPTIONS.FUNC = 13;
OPTIONS.PARALLEL = true;
OPTIONS.REAL = false;
OPTIONS.RADIATION = 0;
OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);
OPTIONS.FOURIER_AUX = [-3*pi/8, 3*pi/8];

run('DRIVER_exact');

errors3 = zeros(5,1);
errors3(1) = errors;

Mc = 40;
M = [Mc, Mf, Mr];
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);
OPTIONS.GRIDS = [7:10];

run('DRIVER_exact');

errors3(2:end) = errors;

ratios3 = log2(errors3(1:end-1) ./ errors3(2:end));


%% output a latex table
n = (2 .^ [6:10])';
output = [n, errors1, [0; (ratios1)], errors2, [0; (ratios2)], errors3, [0; (ratios3)]];
disp(output);

formatspec = "%i & %3.2e & %.2f & %3.2e & %.2f & %3.2e & %.2f \\\\\n";

fprintf(formatspec, output');

