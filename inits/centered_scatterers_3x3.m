% centered_scatterers_3x3.m

CASE = 3;
OPTIONS.Ns = CASE^2;
OPTIONS.types = [2, 2, 2, 2, 1, 2, 2, 2, 2];

SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);

edges = zeros(OPTIONS.Ns,4); start = [-3, -1, -3, -1];
for i=1:CASE
    for j=1:CASE
        edges(3*(i-1)+j,:) = [start(1)+(j-1)*2, start(2)+(j-1)*2, ...
                            start(3)+(i-1)*2, start(4)+(i-1)*2];
    end
end
      
OPTIONS.ADJACENCY = subdomain_adjacency(edges);

k = [3*ones(OPTIONS.Ns,1), 10*ones(OPTIONS.Ns,1)];
k(5,2) = 0; % the one without a scattering rod, center

M = [30*ones(OPTIONS.Ns,1), 35*ones(OPTIONS.Ns,1)];
M(5,2) = 0;