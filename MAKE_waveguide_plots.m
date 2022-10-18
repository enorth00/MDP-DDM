% MAKE_waveguide_plots.m
%
% Makes seven plots. Four plots for the different buffer cases, and three
% error plots.

k = 1.229;
figure();
% ERROR = struct('U', cell(1,OPTIONS.Ns));
% PLOT_VAL=3;
count = 1;
for BUFFS=1:4
    subplot(4,2,count);
    filename = "k_" + k + "_bus_buff" + BUFFS + ".mat";
    load(filename);
    PLOT_waveguide;
    count = count + 1;
    
    if(BUFFS == 1)
        A = SOLUTION;
        count = count + 1;
    else
        B = SOLUTION;
        subplot(4,2,count);
        PLOT_waveguide_errors;
        count = count + 1;
        A = B;
    end
end

% 
% k = 1.229;
% ERROR = struct('U', cell(1,OPTIONS.Ns));
% for BUFFS=1:4
%     if(BUFFS~=1)
%         A = SOLUTION;
%     end
%     filename = "k_" + k + "_bus_buff" + BUFFS + ".mat";
%     load(filename);
%     if(BUFFS~=1)
%         B = SOLUTION;
%         PLOT_waveguide_errors;
%     end
%     PLOT_waveguide;
% end