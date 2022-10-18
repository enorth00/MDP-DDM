% allocate_subdomains.m
% 
% Allocates the reference subdomains
function SUBDOMAIN = allocate_subdomains(Ns)
    % Ns = (N)umber of (s)ubdomains
    SUBDOMAIN = struct('params', cell(1,Ns), ... % k, M, alpha, beta
                       'boundary', cell(1,Ns), ... % functions for BCs
                       'source', cell(1,Ns), ... % function for source
                       'edges', cell(1,Ns), ... % actual coords of edges
                       'type', cell(1,Ns), ... % type to match reference
                       'basis', cell(1,Ns) ... % four chars, c or f
                       );
end