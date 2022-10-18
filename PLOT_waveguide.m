% PLOT_waveguide

% figure(); % comment out if doing subplots

PLOT_VAL = 3;
Ns = OPTIONS.Ns;
width = 15;
buffers = ((Ns/width) - 1) / 2;

SUBS_TO_PLOT = (buffers*width + 1) : (buffers*width + width);

mx=NaN; Mx=NaN; my=NaN; My=NaN; mf = NaN; Mf = NaN;

hold on;
for i = SUBS_TO_PLOT
    
    T = SUBDOMAIN(i).type;
    switch(PLOT_VAL)
        case 1
            uplot = abs(SOLUTION(i).U);
        case 2
            uplot = real(SOLUTION(i).U);
        case 3
            uplot = imag(SOLUTION(i).U);
    end
    
    if(T==2)
        uplot(~(BLOCK(T).grids.Mp | BLOCK(T).grids.rodMp)) = NaN;
    else
        uplot(~(BLOCK(T).grids.Mp)) = NaN;
    end
    
    mx = min(mx, SUBDOMAIN(i).edges(1));
    Mx = max(Mx, SUBDOMAIN(i).edges(2));
    my = min(my, SUBDOMAIN(i).edges(3));
    My = max(My, SUBDOMAIN(i).edges(4));
    mf = min(mf, min(min(uplot)));
    Mf = max(Mf, max(max(uplot)));
    
    xt = @(x) x + (SUBDOMAIN(i).edges(1) - BLOCK(T).edges(1));
    yt = @(y) y + (SUBDOMAIN(i).edges(3) - BLOCK(T).edges(3));

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
% title(title_val + " of Solution, k = " + SUBDOMAIN(1).params.k(1) + ...
%       ", Buffers = " + buffers);
title(buffers + " Buffer Rows");
colorbar;

% xlabel('x');
% ylabel('y');
% zlabel('u');

axis([mx-1, Mx+1, my-1, My+1, mf, Mf]);
caxis([-0.42, 0.42]);
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


