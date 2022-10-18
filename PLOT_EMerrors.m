% PLOT_EMerrors.m

figure();
semilogy(MF,e_vec)
legend('Sommerfeld', 'EM2', 'EM3', 'EM4');
xlabel('M');
ylabel('||Reflection Error||_\infty')
set(gca, 'FontSize', 14)
grid on;
kstr = num2str(k1);
title("Reflection Error of each EM ABC (k=" + kstr + ")");