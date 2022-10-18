% rod_base.m

OPTIONS.T = 2; 

Mc = 30; % number of chebyshev functions per side
Mf = 51; % number of fourier functions per side
Mr = 35; % number of fourier functions per rod

M = [Mc, Mf, Mr];
k1 = 3; k2 = 3;

TYPES = cell(1,3);
TYPES(1,:) = [{1}, {'cccc'}, {[k1, 0]}];
TYPES(2,:) = [{2}, {'cccc'}, {[k1, k2]}];

DOMAIN = [2]; 
start = [-1, 1, -1, 1]; % edges for \Omega_1
%  ___
% | 1 |
% |___|

OPTIONS.Ns = numel(DOMAIN);
OPTIONS.TYPES = TYPES;
BLOCK = allocate_blocks([1, 2]);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);