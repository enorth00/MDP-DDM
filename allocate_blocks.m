% allocate_blocks.m
% 
% Allocates the reference subdomain building blocks
function BLOCK = allocate_blocks(types)
    % types = which reference subdomains are used
    T = max(types);
    BLOCK = struct('edges', cell(1,T), ... % parameters for this reference
                   'grids', cell(1,T), ... % logical grids for this ref
                   'N', cell(1,T), ... % grid dimension for fast access
                   'radius', cell(1,T), ...
                   'aux', cell(1,T), ... % spacing for the auxiliary domain
                   'shape', cell(1,T) ... % identifier for the geometry
                   );
                   
    BLOCK = define_blocks(BLOCK, types);
    
    
end

% This is where the reference parameters for each block are predefined.
function BLOCK = define_blocks(BLOCK, types)
    for i=types
        switch(i)
            case 1
                % standard empty block
                BLOCK(i) = Define_Block1(BLOCK(i));
            case 2
                % standard scattering rod block
                BLOCK(i) = Define_Block2(BLOCK(i));
            case 3
                % Circle to fill in the scattering rod
                BLOCK(i) = Define_Block3(BLOCK(i));
        end
    end

end

function BLOCK = Define_Block1(BLOCK)
    BLOCK.edges = [-1, 1, -1, 1];
    BLOCK.radius = 0;
    BLOCK.aux = 0.2;
    BLOCK.shape = 1;
end

function BLOCK = Define_Block2(BLOCK)
    BLOCK.edges = [-1, 1, -1, 1];
    BLOCK.radius = 0.37;
    BLOCK.aux = 0.2;
    BLOCK.shape = 2;
end

function BLOCK = Define_Block3(BLOCK)
    BLOCK.radius = 0.37;
    BLOCK.aux = 0.2;


end