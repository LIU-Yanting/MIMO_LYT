clc,clear
C_33_double = MIMO_Capacity(2);

C_33_single = MIMO_Capacity(1);

C_33_uncor = MIMO_Capacity(0);
C_33_uncor = ones(length(C_33_double), 1) * C_33_uncor;

C_22_uncor = MIMO_Capacity(0,2,2);
C_22_uncor = ones(length(C_33_double), 1) * C_22_uncor;

C_11_uncor = MIMO_Capacity(0,1,1);
C_11_uncor = ones(length(C_33_double), 1) * C_11_uncor;

ECC = readtable('ECC.xlsx');
%%
figure(1)
plot(ECC.Freq, C_33_double,'-dr', ECC.Freq, C_33_single,'-pg', ECC.Freq, C_33_uncor,'--b', ECC.Freq, C_22_uncor, '-.k',...
    ECC.Freq, C_11_uncor,'-c','LineWidth',2)
title('Comparison of Ergodic Capacity in MIMO and SISO Systems','FontSize',18)
xlabel('Frequency (GHz)','FontSize',18)
ylabel('Ergodic Capacity (bit/s/Hz)','FontSize',18)
set(gca, 'Fontname', 'Times New Roman','FontSize',15);
legend('Tri-polarized Antenna Arrays at Tx and Rx','Tri-polarized Antenna Array at Tx only',...
    'Ideal $3\times3$ MIMO System','Ideal $2\times2$ MIMO System',...
    'Ideal SISO System', 'Interpreter','latex','FontSize',18)
%%
figure(2)
ratio_d = C_33_double./C_11_uncor;
ratio_s = C_33_single./C_11_uncor;
plot(ECC.Freq, ratio_d, '-dr',ECC.Freq, ratio_s,'-pg', 'LineWidth',2)
title('Capacity Enhancement Factors of single- and double-side Tri-polarized Antenna Arrays','FontSize',18)
xlabel('Frequency (GHz)','FontSize',18)
ylabel('Capacity Enhancement Factor','FontSize',18)
set(gca, 'Fontname', 'Times New Roman','FontSize',15);
legend('Tri-polarized Antenna Arrays at Tx and Rx','Tri-polarized Antenna Array at Tx only',...
   'Interpreter','latex','FontSize',18)

