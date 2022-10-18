T = 1;
BLOCK = struct('edges', cell(1,T), ... % parameters for this reference
               'grids', cell(1,T), ... % logical grids for this ref
               'N', cell(1,T), ... % grid dimension for fast access
               'radius', cell(1,T), ...
               'aux', cell(1,T), ... % spacing for the auxiliary domain
               'shape', cell(1,T) ... % identifier for the geometry
               );
BLOCK.edges = [-1, 1, -1, 1];
BLOCK.radius = 0.37;
BLOCK.aux = 0.2;
BLOCK.shape = 2;
BLOCK.N = [17,17];

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

    c = edges; r = rad;
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
    Mm = ~Mp;
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
    
    
    
    
    if(rad>0)
        %     Mm = ~Mp;   % Mm is only needed here to define Mm
    N1 = length(Mp(:,1));  % N1 is length of x, N2 is length of y
    N2 = length(Mp(1,:));
%     Mm(:,1) = false(N1,1);
%     Mm(:,N2) = false(N1,1);
%     Mm(1,:) = false(1,N2);
%     Mm(N1,:) = false(1,N2);
    
    rodNp = false(N1,N2);
    rodNm = false(N1,N2);
    
    for i=2:N1-1
        for j=2:N2-1
            if(rodMp(i,j))
              rodNp(i,j) = true;
              rodNp(i+1,j) = true;
              rodNp(i-1,j) = true;
              rodNp(i,j+1) = true;
              rodNp(i,j-1) = true;
              rodNp(i+1,j+1) = true;
              rodNp(i+1,j-1) = true;
              rodNp(i-1,j+1) = true;
              rodNp(i-1,j-1) = true;
            else
              rodNm(i,j) = true;
              rodNm(i+1,j) = true;
              rodNm(i-1,j) = true;
              rodNm(i,j+1) = true;
              rodNm(i,j-1) = true;
              rodNm(i+1,j+1) = true;
              rodNm(i+1,j-1) = true;
              rodNm(i-1,j+1) = true;
              rodNm(i-1,j-1) = true;
            end
        end
    end
    rodgam = rodNp==rodNm;
    end
