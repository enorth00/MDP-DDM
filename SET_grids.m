function grids = SET_grids(BLOCK)
    % Convenience def for number of gridpoints
    N = BLOCK.N;
    N1 = N(1); 
    N2 = N(2);

    % Establish the x and y grids (for the nodes)
    edges = BLOCK.edges; % typically [-1, 1, -1, 1]
    aux = BLOCK.aux;
    rad = BLOCK.radius;
    
    a = edges(1) - aux;
    b = edges(2) + aux;
    c = edges(3) - aux;
    d = edges(4) + aux;

    x = linspace(a,b,N1); grids.x = x;
    y = linspace(c,d,N2); grids.y = y;

    switch(BLOCK.shape)
        case 1
            grids.Mp = rect_def_Mgrids(edges, grids);
        case 2
            [grids.Mp, grids.rodMp] = scatter_def_Mgrids(edges, rad, grids);
    end
    [grids.Np, grids.Ngamma] = def_Ngrids(grids.Mp);
    if(rad>0)
        [grids.rodNp, grids.rodNgamma] = def_Ngrids(grids.rodMp);
    end
    % Handles "cases" of blocks by use of rad > 0 or not.
    grids.sides = square_sides_assignment(grids.Ngamma, rad);

    [XX1,YY1] = meshgrid(x, y);
    % Meshgrid flips the ordering so that in (i,j) indexing, j is
    % associated with x and i with y, instead of vice-versa like I've been
    % using all along.
    XX = XX1'; YY = YY1';
    grids.coords = [XX(grids.Ngamma), YY(grids.Ngamma), ...
                    grids.sides(grids.Ngamma)];
end


%%%%%%-----Helper Functions-----%%%%%
function Mp = rect_def_Mgrids(c, grids)
    % c is a shortened name for the edges
    x = grids.x; y = grids.y;
    N1 = length(x); N2 = length(y);
    Mp = false(N1, N2);
    for i=1:N1
        for j=1:N2
        % If the x and y coordinates of a node fall
        % *strictly* inside the square, then Mp is true
            if(x(i) > c(1) && x(i) < c(2) && y(j) > c(3) && y(j) < c(4))
                Mp(i,j) = true;
            end
        end
    end
end
function [Mp, rodMp] = scatter_def_Mgrids(c, r, grids)
    % c is a shortened name for the edges
    x = grids.x; y = grids.y;
%     p = [(max(x)+min(x))/2, (max(y)+min(y))/2]; % center of the rod
    N1 = length(x); N2 = length(y);
    Mp = false(N1, N2);
    rodMp = false(N1, N2);
    for i=1:N1
        for j=1:N2
        % If the x and y coordinates of a node fall
        % *strictly* inside the square, then Mp is true
            if(x(i) > c(1) && x(i) < c(2) && y(j) > c(3) && y(j) < c(4))
                Mp(i,j) = true;
            end
            % Remove those nodes that are inside the square.
            if((x(i))^2 + (y(j))^2 <= r^2)
                Mp(i,j) = false;
                rodMp(i,j) = true;
            end
        end
    end
end
function [Np, gam] = def_Ngrids(Mp)
%     Mm = ~Mp;   % Mm is only needed here to define Mm
    N1 = length(Mp(:,1));  % N1 is length of x, N2 is length of y
    N2 = length(Mp(1,:));
%     Mm(:,1) = false(N1,1);
%     Mm(:,N2) = false(N1,1);
%     Mm(1,:) = false(1,N2);
%     Mm(N1,:) = false(1,N2);
    
    Np = false(N1,N2);
    Nm = false(N1,N2);
    
    for i=2:N1-1
        for j=2:N2-1
            if(Mp(i,j))
              Np(i,j) = true;
              Np(i+1,j) = true;
              Np(i-1,j) = true;
              Np(i,j+1) = true;
              Np(i,j-1) = true;
              Np(i+1,j+1) = true;
              Np(i+1,j-1) = true;
              Np(i-1,j+1) = true;
              Np(i-1,j-1) = true;
            else
              Nm(i,j) = true;
              Nm(i+1,j) = true;
              Nm(i-1,j) = true;
              Nm(i,j+1) = true;
              Nm(i,j-1) = true;
              Nm(i+1,j+1) = true;
              Nm(i+1,j-1) = true;
              Nm(i-1,j+1) = true;
              Nm(i-1,j-1) = true;
            end
        end
    end
    gam = Np==Nm;
end
function sides = square_sides_assignment(Ngamma, r)

    [N1, N2] = size(Ngamma);
%     NN = length(Ngamma(Ngamma));
%     M = NN/4;
    sides = zeros(N1,N2);
    counter = 0;
    M = 0; inset = 0;
    % Each side is done individually. The order of each pair
    % of for-loops is different, and we effectively count out
    % 1/4th of the total boundary nodes to be on each side. 
    % This process guarantees an even and consistent distribution
    % of nodes to each side.
    
    for i=1:N1
        if(counter==0)
            for j=1:N2
                if(Ngamma(i,j))
                    % Mark as side 1 and add to the count
                    counter = counter+1;
                    sides(i,j) = 1;

                    % Skip the first time, and don't accidentally
                    % overcount. These conditions help the corners
                    % work out correctly. This action also gives us 
                    % the "two-layers" that each side ends up having.
                    if(counter>1)
                        sides(i+1,j) = 1;
                        counter = counter+1;
                    end
                end
            end
        else
            M = counter - 3;
            break;
        end
    end

    counter = 0;
    for i = N1:-1:1
        for j = N2:-1:1
            if(Ngamma(i,j))
                inset = j;
                if(counter<M)
                    counter = counter+1;
                    sides(i,j) = 2;

                    if(counter>1 && counter < M)
                        sides(i-1,j) = 2;
                        counter = counter+1;
                    end
                end
            end
        end
        if(counter > 0)
            break;
        end
    end

    counter = 0;
    for j = 1:N2
        for i = N1:-1:1
            if(Ngamma(i,j))
                if(counter<M)
                    counter = counter+1;
                    sides(i,j) = 3;

                    if(counter>1 && counter<M)
                        sides(i,j+1) = 3;
                        counter = counter+1;
                    end
                end
            end
        end
        if(counter > 0)
            break;
        end
    end

    counter = 0;
    for j = N2:-1:1
        for i = 1:N1
            if(Ngamma(i,j))
                if(counter < M)
                    counter = counter+1;
                    sides(i,j) = 4;

                    if(counter>1 && counter<M)
                        sides(i,j-1) = 4;
                        counter = counter+1;
                    end
                end
            end
        end
        if(counter>0)
            break;
        end
    end
    
    if(r > 0)
        for i = inset+5:N1-inset-5
            for j = inset+5:N2-inset-5
                if(Ngamma(i,j))
                    sides(i,j) = 5;
                end
            end
        end
    end
end