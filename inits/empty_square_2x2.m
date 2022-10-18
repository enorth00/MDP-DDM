% empty_square_2x2.m

OPTIONS.Ns = 4;
OPTIONS.types = [1, 1, 1, 1];

SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);

edges = [-1, 1, -1, 1;
          1, 3, -1, 1;
         -1, 1,  1, 3;
          1, 3,  1, 3];
      
OPTIONS.ADJACENCY = subdomain_adjacency(edges);


k = 3*[ones(4,1), zeros(4,1)];
M = 31*[ones(4,1), zeros(4,1)];