% tunnel_5x5.m

CASE = 5;
OPTIONS.Ns = CASE^2;
o = ones(1,5);
OPTIONS.types = [2*[o, o], o, 2*[o, o]];

SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);

edges = zeros(OPTIONS.Ns,4); start = [-1, 1, -5, -3];
for i=1:CASE
    for j=1:CASE
        edges(CASE*(i-1)+j,:) = [start(1)+(j-1)*2, start(2)+(j-1)*2, ...
                            start(3)+(i-1)*2, start(4)+(i-1)*2];
    end
end
      
OPTIONS.ADJACENCY = subdomain_adjacency(edges);


k = [3, 10;
     3, 10;
     3, 10;
     3, 10;
     3, 10;
     3, 10;
     3, 10;
     3, 10;
     3, 10;
     3, 10;
     3, 0;
     3, 0;
     3, 0;
     3, 0;
     3, 0;
     3, 10;
     3, 10;
     3, 10;
     3, 10;
     3, 10;
     3, 10;
     3, 10;
     3, 10;
     3, 10;
     3, 10];

% k = [3, 3;
%      3, 3;
%      3, 3;
%      3, 0;
%      3, 0;
%      3, 0;
%      3, 3;
%      3, 3;
%      3, 3];
M = [30, 25;
     30, 25;
     30, 25;
     30, 25;
     30, 25;
     30, 25;
     30, 25;
     30, 25;
     30, 25;
     30, 25;
     30, 0;
     30, 0;
     30, 0;
     30, 0;
     30, 0;
     30, 25;
     30, 25;
     30, 25;
     30, 25;
     30, 25;
     30, 25;
     30, 25;
     30, 25;
     30, 25;
     30, 25];