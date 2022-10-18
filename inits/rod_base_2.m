% rod_base_2.m

OPTIONS.Ns = 2;
OPTIONS.types = [2,2];

SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);

edges = [-1, 1, -1, 1;
          1, 3, -1, 1];
OPTIONS.ADJACENCY = subdomain_adjacency(edges);

k = [3, 10;
     3, 10];
M = [30, 35;
     30, 35];