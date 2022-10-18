% empty_base_2.m

OPTIONS.Ns = 2;
OPTIONS.types = [1,1];

SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);

edges = [-1, 1, -1, 1;
          1, 3, -1, 1];
OPTIONS.ADJACENCY = subdomain_adjacency(edges);

k = [3, 0;
     3, 0];
M = [30, 0;
     30, 0];