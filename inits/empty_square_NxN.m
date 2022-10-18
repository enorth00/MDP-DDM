% empty_square_NxN.m

CASE = 3;
OPTIONS.Ns = CASE^2;
OPTIONS.types = ones(1,OPTIONS.Ns);

SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);

edges = zeros(OPTIONS.Ns,4); start = [-CASE, -CASE+2, -CASE, -CASE+2];
for i=1:CASE
    for j=1:CASE
        edges(CASE*(i-1)+j,:) = [start(1)+(j-1)*2, start(2)+(j-1)*2, ...
                            start(3)+(i-1)*2, start(4)+(i-1)*2];
    end
end
      
OPTIONS.ADJACENCY = subdomain_adjacency(edges);


k = 5*[ones(OPTIONS.Ns,1), zeros(OPTIONS.Ns,1)];
M = [31*ones(OPTIONS.Ns,1), zeros(OPTIONS.Ns,1)];
