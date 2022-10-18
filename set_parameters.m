function SUBDOMAIN = set_parameters(SUBDOMAIN, k, M, edges, bc, type)
% type will determine the nature of the building block
% 0 = 'rectangle'
% 1 = 'rod'

    params.k = k(1);
    params.krod = k(2);
    params.M = M;
    SUBDOMAIN.edges = edges;
    params.alpha = bc(1,:);
    params.beta = bc(2,:);
    SUBDOMAIN.type = type;
    
    SUBDOMAIN.params = params;
end