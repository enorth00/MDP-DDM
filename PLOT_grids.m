DEF_allGrids;
figure();
[XX,YY] = meshgrid(x, y);
h = x(2) - x(1);
% gr = g.rodNp;
% F = gr(gr);
F1 = rodgam;
F2 = Mm;


scatter(XX(F1), YY(F1), 50, 'k', 'filled');
hold on;
% scatter(XX(F2), YY(F2), 50, 'b', 'filled');
% rectangle('Position', [-1, -1, 2, 2], 'FaceColor', 'none', 'EdgeColor', ...
%     'black', 'LineWidth', 1);
circle(0,0,0.37);


num_steps = 10;
xticks(-num_steps*h:h:num_steps*h); yticks(-num_steps*h:h:num_steps*h);
xticklabels({}); yticklabels({});
xlim([-num_steps*h, num_steps*h]); ylim([-num_steps*h, num_steps*h]);
grid on;
set(gca,'XColor','none')
set(gca,'YColor','none')
hold off;