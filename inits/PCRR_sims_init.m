% PCRR_sims_init.m

OPTIONS.T = 2; % this refers to how many reference blocks we will need
                   % to load in ultimately. Each unique k-profile for a
                   % given block shape will require its own.

Mc = 15; % number of chebyshev functions per side
Mf = 21; % number of fourier functions per side. Define this in runner
Mr = 21; % number of fourier functions per rod

M = [Mc, Mf, Mr];

k1 = 1.13; 
% k1 = 1.08;

k2 = k1*3.48;

TYPES = cell(1,3);

TYPES(1,:) = [{1}, {'cccc'}, {[k1, 0]}];
TYPES(2,:) = [{2}, {'cccc'}, {[k1, k2]}];

% RING2 = 2*ones(6); RING2([2,5],2:5) = 1; RING2(2:5,[2,5]) = 1;
RING3 = 2*ones(7); RING3([2,6],2:6) = 1; RING3(2:6,[2,6]) = 1;

% Full size PCRR (25x17, ring3)
%DOMAIN = [2*ones(4,25); ones(1,25); 2*ones(7, 9), RING3, 2*ones(7,9)];
%DOMAIN = [DOMAIN; ones(1,25); 2*ones(4,25)];
%start = [-3, -1, 7, 9];

% Just the bus
height = 1; width = 3;
line = 2*ones(height,width);
DOMAIN = [line; ones(1,width); line; ones(1, width); line];
start = [-3, -1, -1+2*height, 1+2*height];

% Pared-down PCRR (15x13, ring3)
% DOMAIN = [2*ones(2,15); ones(1,15); [2*ones(7,4), RING3, 2*ones(7,4)]; ones(1,15); 2*ones(2,15)];
% start = [-3, -1, 3, 5];

OPTIONS.Ns = numel(DOMAIN);
OPTIONS.TYPES = TYPES;
BLOCK = allocate_blocks([1, 2]);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);
