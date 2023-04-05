clear all 
close all
nmLIB = [400 450 500 550 600 650 700 725];
muawater = [0.0000663,0.0000922,0.000204,0.000565,0.002224,0.0034,0.00624,0.01489];
muadeoxy = [223296 103292 20862 53412 14677.2 3750.12 1794.28 1244.44]; 
muaoxy = [266232 62816	20932.8	43016 3200	368 290	364];
muafat = [15 6.375 1.906 0.773 0.466 0.471 0.323 0.415];
muafat = muafat./100;


%get absorption from extrinction
l = length(nmLIB);
for freq = 1:l
    muadeoxy(freq) = 2.303 * muadeoxy(freq) * 150 / 64500;
    muaoxy(freq) = 2.303 * muaoxy(freq) * 150 / 64500;
end

save spectralLIBv2.mat muawater muafat muadeoxy muaoxy nmLIB 

hold on
plot(nmLIB,muawater,'DisplayName','water');
plot(nmLIB,muaoxy,'DisplayName','hbO2');
plot(nmLIB,muafat,'DisplayName','fat');
plot(nmLIB,muadeoxy,'DisplayName','hb');
legend()
set(gca, 'YScale', 'log')
hold off



load spectralLIB

figure
hold on
plot(nmLIB,muawater,'DisplayName','water');
plot(nmLIB,muaoxy,'DisplayName','hbO2');
plot(nmLIB,muafat,'DisplayName','fat');
plot(nmLIB,muadeoxy,'DisplayName','hb');
legend()
set(gca, 'YScale', 'log')
hold off


nmLIB = [470 530 660 770 810 940 1020 1050];
muawater = [2.47E-05, 3.20E-05, 0.00032, 0.0019858, 0.026737, 0.036, 0.052, 0.0885];
muaoxy = [17.05, 20.91, 0.15, 0.3, 0.4, 0.65, 0.61, 0.49];
muadeoxy = [8.29, 17.82, 1.64, 0.73, 0.45, 0.43, 0.28, 0.22];
muafat = [0.076, 0.073, 0.065, 0.110, 0.138, 0.144, 0.016, 0.157];

muawater = muawater.*10;
muaoxy = muaoxy.*10;
muadeoxy = muadeoxy.*10;
muafat = muafat.*10;



save spectralLIBv3.mat muawater muafat muadeoxy muaoxy nmLIB 

figure
hold on
plot(nmLIB,muawater,'DisplayName','water');
plot(nmLIB,muaoxy,'DisplayName','hbO2');
plot(nmLIB,muafat,'DisplayName','fat');
plot(nmLIB,muadeoxy,'DisplayName','hb');
legend()
set(gca, 'YScale', 'log')
hold off
