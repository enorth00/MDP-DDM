% PLOT_resonator.m

% figure();

PLOT_VAL = 3;
subs_to_plot = [50:56, 65:71, 80:86, 95:101, 110:116, 125:131, 140:146];

mx=NaN; Mx=NaN; my=NaN; My=NaN; mf = NaN; Mf = NaN;

for i = subs_to_plot
    
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
    hold on;
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
title(title_val + " of k = " + num2str(SUBDOMAIN(50).params.k(1)));
colorbar;

xlabel('x');
ylabel('y');
zlabel('u');

axis([mx-1, Mx+1, my-1, My+1, mf, Mf]);
% caxis([-0.42, 0.42]);
set(gca, 'FontSize', 16);

for i=subs_to_plot
    edge = SUBDOMAIN(i).edges;
    rectangle('Position', [edge(1), ...
                           edge(3), ...
                           edge(2)-edge(1), ...
                           edge(4)-edge(3)]);
end
for i=subs_to_plot
    if(BLOCK(SUBDOMAIN(i).type).radius > 0)
        circle((SUBDOMAIN(i).edges(2)+SUBDOMAIN(i).edges(1))/2, ...
               (SUBDOMAIN(i).edges(4)+SUBDOMAIN(i).edges(3))/2, ...
                BLOCK(SUBDOMAIN(i).type).radius);
    end
end

view(0,90);

