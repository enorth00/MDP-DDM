K = [1, 1.03, 1.15, 1.16, 1.294, 1.303];


figure();
for i = 1:length(K)
    load("k_" + num2str(K(i)) + ".mat");
    subplot(3,2,i);
%     PLOT_PCRR;
    PLOT_resonator;
end
