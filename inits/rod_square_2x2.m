% rod_square_2x2.m

OPTIONS.Ns = 4;
OPTIONS.types = [2,2,2,2];

SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);

edges = [-1, 1, -1, 1;
          1, 3, -1, 1;
         -1, 1,  1, 3;
          1, 3,  1, 3];


OPTIONS.ADJACENCY = subdomain_adjacency(edges);

o = ones(OPTIONS.Ns,1);
k = [10*o, 10*o];
M = [30*o, 25*o];