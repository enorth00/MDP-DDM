function [SOLUTION, QR_time, solve_time] = problem_solver(SUBDOMAIN, BLOCK, OPTIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   problem_solver
%
%   Inputs:     
%               
%           
%   Purpose:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = OPTIONS.Ns;
    SOLUTION = struct('c', cell(1,n), 'U', cell(1,n), 'rod', cell(1,n));
    
    % Workhorse. 
%     disp('Source Construction');
    Ft = source_construction(BLOCK, SUBDOMAIN, OPTIONS);

%     disp('Handling boundary conditions');
    SOLUTION = robin_boundaries(SOLUTION, Ft, BLOCK, SUBDOMAIN, OPTIONS);

%     disp('Applying difference potentials');
    SOLUTION = potential_handling(BLOCK, SUBDOMAIN, Ft, SOLUTION, OPTIONS);
    
    % These are here until I decide if I want to record timings.
    solve_times = zeros(1,n);
    QR_time = 0;
    solve_time = sum(solve_times)*0;
end

function Ft = source_construction(BLOCK, SUBDOMAIN, OPTIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   source_construction
%
%   Inputs:     
%               
%           
%   Purpose:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = OPTIONS.Ns;
    Ft = struct('F', cell(1,n), 'GBf', cell(1,n), 'ext', cell(1,n), 'Frod', cell(1,n), 'rodext', cell(1,n));
    
    if(OPTIONS.PARALLEL)
        parfor i=1:n
            T = SUBDOMAIN(i).type;
            [Ft(i).F, Ft(i).GBf, Ft(i).ext, Ft(i).Frod, Ft(i).rodext] = source_extender(BLOCK(T), SUBDOMAIN(i), OPTIONS);
        end
    else
        for i=1:n
            T = SUBDOMAIN(i).type;
            [Ft(i).F, Ft(i).GBf, Ft(i).ext, Ft(i).Frod, Ft(i).rodext] = source_extender(BLOCK(T), SUBDOMAIN(i), OPTIONS);
        end
    end
end
function SOLUTION = robin_boundaries(SOLUTION, Ft, BLOCK, SUBDOMAIN, OPTIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   robin_boundaries
%
%   Inputs:     
%               
%           
%   Purpose:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convenience-defs
    PAR = OPTIONS.PARALLEL;
    adj = OPTIONS.ADJACENCY;
    n = OPTIONS.Ns;
    Q_store = struct('Qstar', cell(1,n), 'Qbar', cell(1,n), 'Qr', cell(1,n));
    
    % load in reference Q-matrices (for each "type")
    N = BLOCK(1).N;
    Q_ref = load_Q_matrices(N(1), SUBDOMAIN, OPTIONS);
    
    % Case out the resolution of boundaries based on the subdomain
    if(PAR)
        parfor i=1:n
            Q_store(i) = boundary_resolution(Q_ref{1,SUBDOMAIN(i).type}, SUBDOMAIN(i), adj(i,:), OPTIONS);
        end
    else
        for i=1:n
            Q_store(i) = boundary_resolution(Q_ref{1,SUBDOMAIN(i).type}, SUBDOMAIN(i), adj(i,:), OPTIONS);
        end
    end

%     disp('Transmission Matching');
%     tic;
    Q = transmission_matching(Q_store, SUBDOMAIN, adj);
%     disp(toc);
    F = source_assembly(Ft, Q_store);
    
%     [Q,F] = add_equations(Q, SUBDOMAIN, F, OPTIONS);
    
%     disp('Beginning QR-factorization');
%     tic;
    c_bar = QR_solution(Q, F);
%     disp(toc);

    SOLUTION = recover_coefficients(c_bar, SOLUTION, SUBDOMAIN, OPTIONS);
end




function SOLUTION = potential_handling(BLOCK, SUBDOMAIN, Ft, SOLUTION, OPTIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   potential_handling
%
%   Inputs:     
%               
%           
%   Purpose:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = OPTIONS.Ns;

    if(OPTIONS.PARALLEL)
        parfor i = 1:n
            % Set the block-index for the relevant reference block. 
            T = SUBDOMAIN(i).type;
            
%             SOLUTION(i).c = DEBUG_get_exact_coeffs(BLOCK(T), SUBDOMAIN(i), OPTIONS);
            % Take coefficients and reconstruct the function handles
            [d0, d2, d4, r] = coeffs_to_functions(SOLUTION(i), SUBDOMAIN(i), OPTIONS);
            
            SOLUTION(i).U = square_potential(d0,d2,d4,r, BLOCK(T), ...
                                             SUBDOMAIN(i), Ft(i), OPTIONS);

            % If we have the rod case, then the rod information needs handled.
            % Extend the circle functions used above, send through N+
            % potential, and truncate to its own M+
            if(T==2)
                SOLUTION(i).rod = rod_potential(r, BLOCK(T), SUBDOMAIN(i), ...
                                                Ft(i), OPTIONS);
%%%%%%%%%%%%%%%%%%%%%%%%% DEBUGGING LOCATION %%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLUTION(i).rod = get_exact(BLOCK(T), SUBDOMAIN(i), BLOCK(T).grids.rodMp, OPTIONS);
% SOLUTION(i).U = get_exact(BLOCK(T), SUBDOMAIN(i), BLOCK(T).grids.Mp, OPTIONS);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                SOLUTION(i).U = SOLUTION(i).U + SOLUTION(i).rod;
            end
            % END ROD CASE
        end
        % END SUBDOMAIN LOOP
    else
        for i = 1:n
            % Set the block-index for the relevant reference block. 
            T = SUBDOMAIN(i).type;
            
%             SOLUTION(i).c = DEBUG_get_exact_coeffs(BLOCK(T), SUBDOMAIN(i), OPTIONS);
            % Take coefficients and reconstruct the function handles
            [d0, d2, d4, r] = coeffs_to_functions(SOLUTION(i), SUBDOMAIN(i), OPTIONS);
            
            SOLUTION(i).U = square_potential(d0,d2,d4,r, BLOCK(T), ...
                                             SUBDOMAIN(i), Ft(i), OPTIONS);

            % If we have the rod case, then the rod information needs handled.
            % Extend the circle functions used above, send through N+
            % potential, and truncate to its own M+
            if(T==2)
                SOLUTION(i).rod = rod_potential(r, BLOCK(T), SUBDOMAIN(i), ...
                                                Ft(i), OPTIONS);
                                            
                %%%%%%%%%%%%%%%%%%%%%%%%% DEBUGGING LOCATION %%%%%%%%%%%%%%%%%%%%%%%%%%
%                 SOLUTION(i).rod = get_exact(BLOCK(T), SUBDOMAIN(i), BLOCK(T).grids.rodMp, OPTIONS);
%                 SOLUTION(i).U = get_exact(BLOCK(T), SUBDOMAIN(i), BLOCK(T).grids.Mp, OPTIONS);    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                SOLUTION(i).U = SOLUTION(i).U + SOLUTION(i).rod;
            end
            % END ROD CASE
        end
        % END SUBDOMAIN LOOP
    end


end

% Called from within robin_boundaries(), these handle the sides of the
% square. Boundaries are resolved with BCs and boundary information, and
% transmission_matching() is the main loop to ultimately 
function Q = boundary_resolution(Q_store, SUBDOMAIN, adj, OPTIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   boundary_resolution
%
%   Inputs:     Qref - 
%          SUBDOMAIN - 
%                adj - 
%
%   Returns:    Q.Qbar - 4 (or 5) cells containing the information of Q
%                        that is "to-be-solved".
%              Q.Qstar - A vector that represents the known information for
%                        set of basis functions that has been multiplied by
%                        the corresponding coefficients. Will be subtracted
%                        from the source contribution.
%               
%   Purpose:    This function is called within a loop over all subdomains.
%               It resolves the boundary conditions, returning "Qstar"
%               which is the vector that compounds the effects of all
%               boundary conditions being sent to the RHS of the BEP, and
%               "Qbar" which is ultimately the columns that account for the
%               unknown boundary data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Needs to return Q.Qstar, Q.Qbar

    % Parameters
    M = SUBDOMAIN.params.M(1:4); Mr = SUBDOMAIN.params.M(5);
    MM = [0, M(1), sum(M(1:2)), sum(M(1:3)), sum(M(1:4))];
%     alpha = SUBDOMAIN.boundary.alpha; beta = SUBDOMAIN.boundary.beta;
    
    [Qm, Qr] = restrict_Q(Q_store, SUBDOMAIN);
    
    % Split the Qref matrix into Q0 and Q1
    Q0 = Qm(:, 1:MM(5)+Mr); Q1 = Qm(:, MM(5)+Mr+1:2*MM(5)+2*Mr);
        
    % Take care of boundary conditions
    [Q_star, Q_bar] = boundary_condition_handler(SUBDOMAIN, Q0, Q1, adj, OPTIONS);

    if(Mr>0)
        Q_bar(5) = {[Q0(:,MM(5)+1:MM(5)+Mr), Q1(:,MM(5)+1:MM(5)+Mr)]};
    end
    
    Q.Qstar = Q_star;
    Q.Qbar = Q_bar;
    Q.Qr = Qr;
end
function Q = transmission_matching(Q_store, SUBDOMAIN, adj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   transmission_matching
%
%   Inputs:     Q_store - Contains the Qbar{1:4} (5 with a rod) that are
%                         used to define the "unknown" parts of the big
%                         matrix Q. These columns correspond to the
%                         coefficients that are getting solved for.
%             SUBDOMAIN - Carries the number of basis functions in
%                         SUBDOMAIN.params.M. We need the whole SUBDOMAIN
%                         structure though because certain SUBDOMAINs need
%                         passed to update_Q_trans
%                   adj - Conveniently passed in separately for the sake
%                         of checking immediately inside the loops.
%           
%   Purpose:    Runs after boundary_resolution. This function handles
%               implanting the columns corresponding to interface sides.
%               For each interface, it checks if its matching side has been
%               accounted for yet. If it has been, then the corresponding
%               columns are added underneath the existing columns. If the
%               matching side has not been accounted for, then this side's
%               columns are simply added as the next section (diagonally).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function will go through each subdomain and insert its relevant
    % information into the "master matrix" Q (which will ultimately be
    % QR'd). For each side of each subdomain, we check to see if the side
    % is exterior or interior, and add columns to the appropriate
    % locations.
    Q = [];
    
    
    n = length(SUBDOMAIN);
    Mskip = zeros(n,5);
    for i=1:n
        M = SUBDOMAIN(i).params.M(1:4);
        Mr = SUBDOMAIN(i).params.M(5);
        % Q_in is the new columns being introduced by SUBDOMAIN(i). These
        % are either the remaining Q_bar columns from an exterior side, or
        % all of the columns (Q0 and Q1) from an interior side whose
        % neighbor is not yet present.
        Q_in = [];
        
        % Each subdomain step should update Q as follows:
        % Q = [Q      , 0   ]
        %     [Q_trans, Q_in]
        
        [n1,n2] = size(Q);
        [m1,~] = size(Q_store(i).Qbar{1});
        Q_trans = zeros(m1,n2);
        for j = 1:4
            if(adj(i,j) == 0)
                Q_in = [Q_in, Q_store(i).Qbar{j}];
                [~, Mskip(i,j)] = size(Q_store(i).Qbar{j});
%                 Mskip(i,j) = M(j);
            else
                if(adj(i,j) > i)
                    % The neighboring subdomain has not been accounted for
                    % yet. All information needs included here.
                    Q_in = [Q_in, Q_store(i).Qbar{j}];
                    Mskip(i,j) = 2*M(j);
                else
                    Q_trans = update_Q_trans(Q_trans, ...
                                             SUBDOMAIN(adj(i,j)), j, ...
                                             Q_store(i).Qbar{j}, ...
                                             Mskip(1:adj(i,j),:));
                end
            end
        end

        
        if(Mr>0) % account for a scattering rod
            Q_in = [Q_in, Q_store(i).Qbar{5}];
            [~,Ms] = size(Q_store(i).Qbar{5});
            Mskip(i,5) = Ms;
        end
        
        if(isempty(Q_in))
            Q = [Q; Q_trans];
        else
            [~,m2] = size(Q_in);
            Q = [Q, zeros(n1,m2); Q_trans, Q_in];
        end
        
        if(Mr>0)
            [a1,b1] = size(Q_store(i).Qr);
            [~,b2] = size(Q);
            Q = [Q; zeros(a1, b2-b1), Q_store(i).Qr];
        end
    end
    
end

% Assembles the (wide) block that will contain the information
% corresponding to those sides that "match up with" existing sides.
function Q_trans = update_Q_trans(Q_trans,SUBDOMAIN,j,Qbar,Mskip)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   update_Q_trans
%
%   Inputs:     
%               
%           
%   Purpose:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MM = sum(sum(Mskip(1:end-1,:)));

    % Loop over "potential targets" within the neighboring subdomains. If
    % its not the side we are neighboring with the current selection, add
    % the corresponding columns to MM and move to the next one.
    ID = [2, 1, 4, 3];
    for z = 1:4
        M = SUBDOMAIN.params.M(z);
        if(ID(z) == j)
            Q_trans(:,MM+1:MM+2*M) = [Qbar(:,1:M), -Qbar(:,M+1:2*M)];
            return;
        end
        MM = MM + Mskip(end,z);
    end

end


%%%--------------- Single-use Stub functions ---------------%%%
% 
% These are functions that are pretty much called once. Logically, it would
% make sense to not separate these into functions and instead call them
% from directly inside the main loop, but they are collectively just lines
% of code that are each trying to accomplish one task. This organization
% keeps this code organized and out of the way during debugging and fixing.
%
%%%---------------------------------------------------------%%%
function Q_ref = load_Q_matrices(N, SUBDOMAIN, OPTIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   load_Q_matrices
%
%   Inputs:     N - The number of grid nodes in each dimension. Provides
%                   the number needed for the name of Q
%       SUBDOMAIN - 
%         OPTIONS - Carries OPTIONS.TYPES and OPTIONS.T. T is the number of
%                   reference blocks being loaded in; one for each unique
%                   geometry/wavenumber profile pairing. TYPES is the
%                   structure defined at init containing the "shape" (1 or
%                   2), and the wavenumber profile in the first and third
%                   spots.
%           
%   Purpose:    Loads in reference Q-matrices. Takes some juggling because
%               we need to make sure all the necessary matrices get 
%               accounted for without loading in superfluous or repetitive
%               information.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OPTIONS.TYPES(*,1) = "empty" or "rod"
    % OPTIONS.TYPES(*,3) = [k, krod]
    T = OPTIONS.TYPES;
    Q_ref = cell(2,OPTIONS.T);
    
    for i=1:OPTIONS.T
        if(T{i,1} == 1)
            shape = "empty";
        elseif(T{i,1} == 2)
            shape = "rod";
        end
        filename = "Q_" + num2str(T{i,3}(1)) + ...
                    "_" + num2str(T{i,3}(2)) + "_" + num2str(N-1);
        Q_ref{1,i} = load("Q_mats/" + shape + "/chebyshev/" + filename + ".mat", 'Q_store');
%         Q_ref{2,i} = load("Q_mats/" + shape + "/fourier/" + filename + ".mat");
    end
    
end
function [Qm, Qr] = restrict_Q(Qin, SUBDOMAIN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   restrict_Q
%
%   Inputs:     Qin - 2x1 cell that contains the reference Q-matrices for
%                     the Chebyshev and Fourier basis functions in the {1}
%                     cell and {2} cell, respectively.
%         SUBDOMAIN - Needed for SUBDOMAIN.params.M and SUBDOMAIN.basis.
%                     params.M should be included as [M1, M2, M3, M4, Mr]
%                     and SUBDOMAIN.basis should be 'ccfc' or similar, 4
%                     characters using 'c' and 'f'.
%           
%   Purpose:    A stub function that is called in boundary_resolution.m.
%               This function accepts the full reference Q-matrices and
%               pares them down to only include the columns needed for our
%               problem.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Q = Qin.Q_store.Q;
    M = SUBDOMAIN.params.M(1:4); Mr = SUBDOMAIN.params.M(5);
    
    % Calculate MMax and MrMax directly from the stored data. This heavily
    % relies on Chebyshev and Fourier reference matrices using the same
    % MMax.
    if(Mr>0)
        [~,Mrtemp] = size(Qin.Q_store.Qrod);
    else
        Mrtemp = 0;
    end
    [NN, Mtemp] = size(Q);
    MrMax = Mrtemp/2;
    MMax = ((Mtemp/2) - MrMax) / 4;
    
    % Allocate a place for the restricted columns to be stored (and
    % returned).
%     [NN, ~] = size(Q{1}.Q);
    Qm = zeros(NN, 2*sum(M));
    
    % MM serves as a "skipping point" so each side knows which column to
    % start at for its own information.
    MM = [0, M(1), sum(M(1:2)), sum(M(1:3)), sum(M)];
    
    % Run the loop over the four sides. Split the restriction formulas
    % based on if the side is chebyshev or fourier.
    for i=1:4
        
        % Chebyshev is straightforward. Take the first M(i) from each
        % block, the rest is ignored.
        if(SUBDOMAIN.basis(i) == 'c')
            Qm(:, MM(i)+1:MM(i+1)) = Q(:, (i-1)*MMax+1 : (i-1)*MMax+M(i));
            Qm(:, MM(5)+Mr+MM(i)+1 : MM(5)+Mr+MM(i+1)) = ...
                Q(:, 4*MMax+MrMax + (i-1)*MMax+1 : 4*MMax+MrMax + (i-1)*MMax + M(i));
        
        % For fourier, the "wasted" information is the tails of the block,
        % we want stuff in the middle. Use a and b as shorthand.
        elseif(SUBDOMAIN.basis(i) == 'f')
            a = (MMax - M(i) + 2)/2; b = (MMax + M(i))/2;
            Qm(:, MM(i)+1:MM(i+1)) = Q{2}.Q(:, (i-1)*MMax + a: (i-1)*MMax + b);
            Qm(:, MM(5)+MM(i)+1:MM(5)+MM(i+1)) = ...
               Q{2}.Q(:, 4*MMax+MrMax + (i-1)*MMax+a : 4*MMax+MrMax + (i-1)*MMax+b);
        end
        
    end
    
    % Included for when there are columns for a scatterer. Repeat of the
    % strategy for fourier functions.
    if(Mr>0)
        a = (MrMax - Mr + 2)/2;
        b = (MrMax + Mr)/2;
        
        % This stuff gets tacked onto the big Q for info in the square
        Qm(:, MM(5)+1:MM(5)+Mr) = Q(:, 4*MMax+a:4*MMax+b);
        Qm(:,2*MM(5)+Mr+1:2*MM(5)+2*Mr) = Q(:,8*MMax+MrMax+a:8*MMax+MrMax+b);
        
        % Internal rod information
        Qr(:,1:Mr) = Qin.Q_store.Qrod(:, a:b);
        Qr(:,Mr+1:2*Mr) = Qin.Q_store.Qrod(:, MrMax+a:MrMax+b);
    else
        Qr = [];
    end
end
function F = source_assembly(Ft, Q_store)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   source_assembly
%
%   Inputs:     Qin - 2x1 cell that contains the reference Q-matrices for
%               
%           
%   Purpose:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NN = 0;
    n = length(Ft);
    
    % NN will be the total length of the RHS vector F.
    for i = 1:n
        NN = NN + length(Ft(i).F);
        if(~isempty(Q_store(i).Qr))
            [n1,~] = size(Q_store(i).Qr);
            NN = NN + n1;
        end
        
    end
    F = zeros(NN,1);
    
    % Nc is a counter for NN of sorts.
    Nc = 0;
    for i=1:n
        NN = length(Ft(i).F);
        F(Nc+1:Nc+NN) = Ft(i).F - Q_store(i).Qstar;
        Nc = Nc + NN;
        if(~isempty(Q_store(i).Qr))
            [n1,~] = size(Q_store(i).Qr);
            F(Nc+1:Nc+n1) = Ft(i).Frod;
            Nc = Nc + n1;
        end
    end
end
function c = QR_solution(Q, F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   QR_solution
%
%   Inputs:     Q - lefthand side of the completed BEP
%               F - righthand side of the completed BEP
%           
%   Purpose:    This serves as a separate housing for computing the QR
%               factorization of the system. This provides a convenient,
%               isolated location in the code where the internal solver can
%               be swapped out easily.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [q,r] = qr(Q, 0);
    c = (r\q')*F;
end
function [d0, d2, d4, r] = coeffs_to_functions(SOLUTION, SUBDOMAIN, OPTIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   coeffs_to_functions
%
%   Inputs:     SOLUTION - contains the solution coefficients, c
%              SUBDOMAIN - contains M (number of basis functions for sides)
%                OPTIONS - carries the interval that external fourier
%                          functions are being expanded over.
%               
%   Purpose:    Given the solution coefficients, this function returns
%               function handles for the evaluation of the boundary
%               functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Convert coefficients into function handles
    c = SOLUTION.c;
    
    % coefficients are stored inside SUBDOMAIN as separate columns
    c_0 = c(:,1);
    c_1 = c(:,2);
    
    % basis functions for the square vs. the rod
    M = SUBDOMAIN.params.M(1:4);
    MM = [0, M(1), sum(M(1:2)), sum(M(1:3)), sum(M(1:4))];
    Mr = SUBDOMAIN.params.M(5);
    
    xi0 = cell(1,4);
    f_basis = OPTIONS.FOURIER_AUX;
    % Allocation of the returned variables.
    d0.xi0 = xi0; d0.xi1 = xi0;
    d2.xi0 = xi0; d2.xi1 = xi0;
    d4.xi0 = xi0;
    for i=1:4
        if(SUBDOMAIN.basis(i) == 'c')
            d0.xi0{i} = chebfun(c_0(MM(i)+1:MM(i+1)), [-1, 1], 'coeffs');
            d0.xi1{i} = chebfun(c_1(MM(i)+1:MM(i+1)), [-1, 1], 'coeffs');
        elseif(SUBDOMAIN.basis(i) == 'f')
            d0.xi0{i} = chebfun(c_0(MM(i)+1:MM(i+1)), f_basis, 'trig', 'coeffs');
            d0.xi1{i} = chebfun(c_1(MM(i)+1:MM(i+1)), f_basis, 'trig', 'coeffs');
        end
        
        d2.xi0{i} = diff(d0.xi0{i}, 2);
        d2.xi1{i} = diff(d0.xi1{i}, 2);

        d4.xi0{i} = diff(d2.xi0{i}, 2);
    end
    if(Mr>0)
        r.d0.xi0{1} = chebfun(c_0(MM(5)+1:MM(5)+Mr), [-pi, pi], 'trig', 'coeffs');
        r.d0.xi1{1} = chebfun(c_1(MM(5)+1:MM(5)+Mr), [-pi, pi], 'trig', 'coeffs');
        r.d2.xi0{1} = diff(r.d0.xi0{1}, 2);
        r.d2.xi1{1} = diff(r.d0.xi1{1}, 2);
        r.d4.xi0{1} = diff(r.d2.xi0{1}, 2);
    else
        r = [];
    end
end
function SOLUTION = recover_coefficients(c_bar, SOLUTION, SUBDOMAIN, OPTIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function:   recover_coefficients
%
%   Inputs:     c_bar - The coefficients solved for in the BEP
%            SOLUTION - Contains space allocated for all coefficients for
%                       all of the subdomains.
%           SUBDOMAIN - Carries M, and gets transferred to
%                       recover_coefficient_handler.
%               
%   Purpose:    Given the coefficients that were solved for in the complete
%               linear system, we can recover the other half of the
%               coefficients using their original substitutions. In the
%               case of external sides, we distribute to the function
%               recover_coefficient_handler to differentiate between
%               interior problems and various external ABCs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    adj = OPTIONS.ADJACENCY;
    [n, ~] = size(adj);
    ID = [2, 1, 4, 3]; % side j is always adjacent to side ID(j)
    for i=1:n
        M = SUBDOMAIN(i).params.M(1:4);
        MM = [0, M(1), sum(M(1:2)), sum(M(1:3)), sum(M)];
        Mr = SUBDOMAIN(i).params.M(5);
        
        
        
        for j=1:4
            if(adj(i,j) == 0)
                SOLUTION(i).c(MM(j)+1:MM(j+1),1:2) = ...
                    recover_coefficient_handler(SUBDOMAIN(i), j, c_bar(1:M(j)), OPTIONS);
                c_bar(1:M(j)) = [];
            else
                % account for interface case
                if(adj(i,j) > i)
                    % This is the case where the info is actually stored
                    SOLUTION(i).c(MM(j)+1:MM(j+1),1:2) = ...
                                            [c_bar(1:M(j)), c_bar(M(j)+1:2*M(j))];
                    c_bar(1:2*M(j)) = [];
                else
                    % The stored info already got processed into an earlier
                    % subdomain. We must retrieve that info from said
                    % subdomain.
                    M2 = SUBDOMAIN(adj(i,j)).params.M(1:4);
                    MM2 = [0, M2(1), sum(M2(1:2)), sum(M2(1:3)), sum(M2(1:4))];
                    
                    SOLUTION(i).c(MM(j)+1:MM(j+1), 1:2) = ...
                         [SOLUTION(adj(i,j)).c(MM2(ID(j))+1:MM2(ID(j)+1), 1), ...
                         -SOLUTION(adj(i,j)).c(MM2(ID(j))+1:MM2(ID(j)+1), 2)];
                end
            end
               
        end
        % If Mr>0, then we strip coefficients into their placeholders.
        % otherwise, do nothing because no coefficients will be present for
        % the rod section.
        if(Mr>0)
            SOLUTION(i).c(MM(5)+1:MM(5)+Mr, 1) = c_bar(1:Mr);
            SOLUTION(i).c(MM(5)+1:MM(5)+Mr, 2) = c_bar(Mr+1:2*Mr);
            
            c_bar(1:2*Mr) = [];
        end
    end

end
function rod = rod_potential(r, BLOCK, SUBDOMAIN, Ft, OPTIONS)

    xi_gamma.Rod = extend_rod_basisfuncs(r, BLOCK, SUBDOMAIN, 0) + ...
                                  extend_rod_basisfuncs(r, BLOCK, SUBDOMAIN, 1) + ...
                                  Ft.rodext;
                
    xi_gamma.XIrod = zeros(BLOCK.N(1), BLOCK.N(2));
    xi_gamma.XIrod(BLOCK.grids.rodNgamma) = xi_gamma.Rod;

    %%%%%%%%%%%%%%%%%%%%%%%%% DEBUGGING LOCATION %%%%%%%%%%%%%%%%%%%%%%%%%%
%     xi_gamma.XIrod = get_exact(BLOCK,SUBDOMAIN,BLOCK.grids.rodNgamma,OPTIONS);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rod = N_plus_potential_rod(xi_gamma.XIrod, BLOCK, SUBDOMAIN, ...
                                           OPTIONS.REAL);
    % Add GBf from within the rod to the potential solution
    rodGBf = get_GBf(BLOCK, SUBDOMAIN, OPTIONS.REAL, true);

    rod = rod + rodGBf;
    rod(~BLOCK.grids.rodMp) = 0;

end
function [U, solve_times] = square_potential(d0,d2,d4,r, BLOCK, SUBDOMAIN, Ft, OPTIONS)
% This returns a VECTOR of values along Ngamma
    xi = extend_xi(d0, d2, d4, r, BLOCK, SUBDOMAIN, Ft.ext);

    % Allocate space for an N1xN2 matrix, then populate those nodes
    %   along Ngamma with values from xi_gamma(i).xi
    XI = zeros(BLOCK.N(1), BLOCK.N(2));
    XI(BLOCK.grids.Ngamma) = xi;

    %%%%%%%%%%%%%%%%%%%%%%%%% DEBUGGING LOCATION %%%%%%%%%%%%%%%%%%%%%%%%%%
%     XI = get_exact(BLOCK, SUBDOMAIN, BLOCK.grids.Ngamma, OPTIONS);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Apply the potential operator P_N+ and add the source term
    [U, solve_times] = N_plus_potential(XI, BLOCK, SUBDOMAIN, OPTIONS.REAL);
    U = U + Ft.GBf;

    % Truncate to M+ (this either gets returned or combined with the
    % information from the rod.
    U(~BLOCK.grids.Mp) = 0;
end
