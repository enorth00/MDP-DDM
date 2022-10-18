% test_domain.m

OPTIONS.T = 1; % this refers to how many reference blocks we will need
                   % to load in ultimately. Each unique k-profile for a
                   % given block shape will require its own.

Mc = 30; % number of chebyshev functions per side
Mf = 31; % number of fourier functions per side
Mr = 35; % number of fourier functions per rod

M = [Mc, Mf, Mr];
k1 = 16; k2 = 0;
%%% Leave room to define each "building block" for the sake of the domain
%%% construction. [{block type}, {edge basis functions}, {wavenumbers}}
TYPES = cell(1,3);
TYPES(1,:) = [{1}, {'cccc'}, {[k1, 0]}];

DOMAIN = [1]; 
start = [-1, 1, -1, 1]; % edges for \Omega_1
%  ___
% | 1 |
% |___|


OPTIONS.Ns = numel(DOMAIN);
OPTIONS.TYPES = TYPES;
BLOCK = allocate_blocks([1, 2]);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);