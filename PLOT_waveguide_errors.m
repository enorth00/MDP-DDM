% PLOT_waveguide_errors.m

analyze_bus_errors;

% figure(); % comment out if doing subplots

PLOT_VAL = 3;

SUBS_TO_PLOT = (buffers*width + 1) : (buffers*width + width);

mx=NaN; Mx=NaN; my=NaN; My=NaN; mf = NaN; Mf = NaN;

hold on;
for i = 1:15
    j = SUBS_TO_PLOT(i);
    T = SUBDOMAIN(i).type;
    switch(PLOT_VAL)
        case 1
            uplot = abs(ERROR(i).U);
        case 2
            uplot = real(ERROR(i).U);
        case 3
            uplot = abs(imag(ERROR(i).U));
    end
    
    if(T==2)
        uplot(~(BLOCK(T).grids.Mp | BLOCK(T).grids.rodMp)) = NaN;
    else
        uplot(~(BLOCK(T).grids.Mp)) = NaN;
    end
    
    mx = min(mx, SUBDOMAIN(j).edges(1));
    Mx = max(Mx, SUBDOMAIN(j).edges(2));
    my = min(my, SUBDOMAIN(j).edges(3));
    My = max(My, SUBDOMAIN(j).edges(4));
    mf = min(mf, min(min(uplot)));
    Mf = max(Mf, max(max(uplot)));
    
    xt = @(x) x + (SUBDOMAIN(j).edges(1) - BLOCK(T).edges(1));
    yt = @(y) y + (SUBDOMAIN(j).edges(3) - BLOCK(T).edges(3));

    [XX, YY] = meshgrid(xt(BLOCK(T).grids.x), yt(BLOCK(T).grids.y));
    
    mesh(XX', YY', uplot);
end

grid on;
switch(PLOT_VAL)
    case 1
        title_val = "Absolute Value";
    case 2
        title_val = "Real Part";
    case 3
        title_val = "Imaginary Part";
end
% title(title_val + " of Error, k = " + SUBDOMAIN(1).params.k(1) + ...
%       ", Buffers = " + (BUFFS-1) + " to " + BUFFS);
title("Error between " + (BUFFS-1) + " and " + (BUFFS) + " Buffers");
colorbar;

% xlabel('x');
% ylabel('y');
% zlabel('u');

axis([mx-1, Mx+1, my-1, My+1, mf, Mf]);
caxis([0, 0.07]);
set(gca, 'FontSize', 16);


% for i=SUBS_TO_PLOT
%     edge = SUBDOMAIN(i).edges;
%     rectangle('Position', [edge(1), ...
%                            edge(3), ...
%                            edge(2)-edge(1), ...
%                            edge(4)-edge(3)]);
% end
% for i=1:Ns
%     if(BLOCK(SUBDOMAIN(i).type).radius > 0)
%         circle((SUBDOMAIN(i).edges(2)+SUBDOMAIN(i).edges(1))/2, ...
%                (SUBDOMAIN(i).edges(4)+SUBDOMAIN(i).edges(3))/2, ...
%                 BLOCK(SUBDOMAIN(i).type).radius);
%     end
% end

hold off;


