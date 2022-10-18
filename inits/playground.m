% test_domain.m

OPTIONS.T = 2; % this refers to how many reference blocks we will need
                   % to load in ultimately. Each unique k-profile for a
                   % given block shape will require its own.

Mc = 25; % number of chebyshev functions per side
Mf = 35; % number of fourier functions per side
Mr = 31; % number of fourier functions per rod

M = [Mc, Mf, Mr];
k1 = 5; k2 = 5;
%%% Leave room to define each "building block" for the sake of the domain
%%% construction. [{block type}, {edge basis functions}, {wavenumbers}}
TYPES = cell(1,3);

TYPES(1,:) = [{1}, {'cccc'}, {[k1, 0]}];

TYPES(2,:) = [{2}, {'cccc'}, {[k1, k2]}];
% TYPES(2,:) = [{1}, {'fccc'}, {[k1, 0]}];

% TYPES(1,:) = [{1}, {'fccf'}, {[k1,0]}];
% TYPES(2,:) = [{1}, {'cfcf'}, {[k1,0]}];
% TYPES(3,:) = [{1}, {'fcfc'}, {[k1,0]}];
% TYPES(4,:) = [{1}, {'cffc'}, {[k1,0]}];

% TYPES(1,:) = [{1}, {'cffc'}, {[k1,0]}];
% TYPES(2,:) = [{1}, {'fffc'}, {[k1,0]}];
% TYPES(3,:) = [{1}, {'fcfc'}, {[k1,0]}];
% TYPES(4,:) = [{1}, {'cfff'}, {[k1,0]}];
% TYPES(5,:) = [{1}, {'ffff'}, {[k1,0]}];
% TYPES(6,:) = [{1}, {'fcff'}, {[k1,0]}];
% TYPES(7,:) = [{1}, {'cfcf'}, {[k1,0]}];
% TYPES(8,:) = [{1}, {'ffcf'}, {[k1,0]}];
% TYPES(9,:) = [{1}, {'fccf'}, {[k1,0]}];

% Fourier Exterior - Chebyshev interior
% TYPES(1,:) = [{1}, {'fccf'}, {[k1,0]}];
% TYPES(2,:) = [{1}, {'cccf'}, {[k1,0]}];
% TYPES(3,:) = [{1}, {'cfcf'}, {[k1,0]}];
% TYPES(4,:) = [{1}, {'fccc'}, {[k1,0]}];
% TYPES(5,:) = [{1}, {'cccc'}, {[k1,0]}];
% TYPES(6,:) = [{1}, {'cfcc'}, {[k1,0]}];
% TYPES(7,:) = [{1}, {'fcfc'}, {[k1,0]}];
% TYPES(8,:) = [{1}, {'ccfc'}, {[k1,0]}];
% TYPES(9,:) = [{1}, {'cffc'}, {[k1,0]}];


%%% Use the previously defined "TYPES" blocks to define a domain in an
%%% array (the ordering will always start with \Omega_1 as the first
%%% element.
DOMAIN = [1,2];
start = [-1, 1, -1, 1];
% start = [-3, -1, 1, 3];
% start = [-5, -3, 3, 5];
%  ___
% | 1 |
% |___|
% | 2 |
% |___|

OPTIONS.Ns = numel(DOMAIN);
OPTIONS.TYPES = TYPES;
BLOCK = allocate_blocks([1, 2]);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);
