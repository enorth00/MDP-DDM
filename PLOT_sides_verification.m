close all;
figure(1);
scatter(0,0,1,0);
hold on;

c = ['r', 'b', 'y', 'm', 'k'];

for i=1:OPTIONS.Ns
    base = BLOCK(SUBDOMAIN(i).type).edges;
    trans = SUBDOMAIN(i).edges;
    
    xt = @(x) x + (trans(1) - base(1));
    yt = @(y) y + (trans(3) - base(3));
    
    coords = BLOCK(SUBDOMAIN(i).type).grids.coords;
    x = xt(coords(:,1)); y = yt(coords(:,2)); z = coords(:,3);
    for j=1:5
        scatter(x(z==j), y(z==j), 16, c(j)); 
    end
end
hold off;
grid on;
legend('0','1', '2', '3', '4', '5');


% scatter(coord(:,1), coord(:,2), 16, coord(:,3), 'filled');



% scatter(x(z==1), y(z==1), 16, c(1), 'filled');
% scatter(x(z==2), y(z==2), 16, c(2), 'filled');
% scatter(x(z==3), y(z==3), 16, c(3), 'filled');
% scatter(x(z==4), y(z==4), 16, c(4), 'filled');
% scatter(x(z==5), y(z==5), 16, c(5), 'filled');
% hold off;
% grid on;
