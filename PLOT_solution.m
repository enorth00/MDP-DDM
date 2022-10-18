% PLOT_solution.m

% close all;
Ns = OPTIONS.Ns;

figure();
hold on;
mx=NaN; Mx=NaN; my=NaN; My=NaN; mf = NaN; Mf = NaN;

for i=1:Ns
    T = SUBDOMAIN(i).type;
%     uplot = imag(SOLUTION(i).U);
%     uplot = abs(error{i});
    uplot = real(SOLUTION(i).U);
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
title('Total Field');
xlabel('x');
ylabel('y');
zlabel('u');
if(mf == 0)
    mf = -1;
end
if(Mf == 0)
    Mf = 1;
end
axis([mx-1, Mx+1, my-1, My+1, mf, Mf]);

for i=1:Ns
    edge = SUBDOMAIN(i).edges;
    rectangle('Position', [edge(1), ...
                           edge(3), ...
                           edge(2)-edge(1), ...
                           edge(4)-edge(3)]);
end
for i=1:Ns
    if(BLOCK(SUBDOMAIN(i).type).radius > 0)
        circle((SUBDOMAIN(i).edges(2)+SUBDOMAIN(i).edges(1))/2, ...
               (SUBDOMAIN(i).edges(4)+SUBDOMAIN(i).edges(3))/2, ...
                BLOCK(SUBDOMAIN(i).type).radius);
    end
end

hold off;




