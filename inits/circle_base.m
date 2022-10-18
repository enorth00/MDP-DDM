% circle_base.m

OPTIONS.Ns = 1;
OPTIONS.types = [3];

SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);

edges = [-1, 1, -1, 1];

OPTIONS.ADJACENCY = subdomain_adjacency(edges);

k = [3, 10];
M = [40, 25];