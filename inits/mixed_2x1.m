% mixed_2x1.m

% How many subdomains are there?
OPTIONS.Ns = 9;
% OPTIONS.types = [2, 2, 2, 1, 1, 1, 2, 2, 2];

% How many different kinds of blocks do we need?
T = 2;
B = cell(4,2);

% What kinds of blocks do I want defined?
B(:,1) = {"empty", 'ccfc', [3, 0], [30, 30, 31, 30, 0]};
B(:,2) = {"rod", 'cccf', [3, 10], [30, 30, 30, 31, 35]};

% Using these sequentially defined blocks, arrange them into a domain
DOMAIN = [1;2];


BLOCK = allocate_blocks([1,2]);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);

edges = zeros(OPTIONS.Ns,4); start = [-1, 1, -3, -1];
for i=1:CASE
    for j=1:CASE
        edges(3*(i-1)+j,:) = [start(1)+(j-1)*2, start(2)+(j-1)*2, ...
                            start(3)+(i-1)*2, start(4)+(i-1)*2];
    end
end
      
OPTIONS.ADJACENCY = subdomain_adjacency(edges);

