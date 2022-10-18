% test_domain.m

OPTIONS.T = 2; % this refers to how many reference blocks we will need
                   % to load in ultimately. Each unique k-profile for a
                   % given block shape will require its own.

% Mc = 40; % number of chebyshev functions per side
Mf = 21; % number of fourier functions per side. Define this in runner
Mr = 21; % number of fourier functions per rod

M = [Mc, Mf, Mr];
k1 = 1.3; k2 = k1*3.48;

TYPES = cell(1,3);

TYPES(1,:) = [{1}, {'cccc'}, {[k1, 0]}];
TYPES(2,:) = [{2}, {'cccc'}, {[k1, k2]}];

DOMAIN = 2*ones(3); DOMAIN(2,2) = 1;
start = [-3, -1, 1, 3];

% DOMAIN = [1];
% start = [-1, 1, -1, 1];

% DOMAIN = [1, 2];
% start = [-1, 1, -1, 1];

% DOMAIN = [2*ones(1,3); ones(1,3); 2*ones(1,3)];
% start = [-3, -1, 1, 3];

%TYPES(1,:) = [{1}, {'fccf'}, {[k1,0]}];
%TYPES(2,:) = [{1}, {'cccf'}, {[k1,0]}];
%TYPES(3,:) = [{1}, {'cfcf'}, {[k1,0]}];
%TYPES(4,:) = [{1}, {'fccc'}, {[k1,0]}];
%TYPES(5,:) = [{1}, {'cccc'}, {[k1,0]}];
%TYPES(6,:) = [{1}, {'cfcc'}, {[k1,0]}];
%TYPES(7,:) = [{1}, {'fcfc'}, {[k1,0]}];
%TYPES(8,:) = [{1}, {'ccfc'}, {[k1,0]}];
%TYPES(9,:) = [{1}, {'cffc'}, {[k1,0]}];

%DOMAIN = [1,2,3; 4,5,6; 7,8,9];
%start = [-3, -1, 1, 3];

OPTIONS.Ns = numel(DOMAIN);
OPTIONS.TYPES = TYPES;
BLOCK = allocate_blocks([1, 2]);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);
