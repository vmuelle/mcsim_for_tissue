close all;
clear all;



%%%%
%This file visualizes the penetretion-depth (PD) and depth-origin (DO) of
%light in skin tissue. The outputfiles of MCML are used for this. 
%%%%

PLOTON = 1;%If PLOTON = 1, the Flux curves are shown for each input file. If PLOTON = 0 only the final plot over all Wavelength is shown.
PD = zeros(1,8); % Parameter for penetration-depth
DO = zeros(1,8); % Parameter for depth-origin 
lambda = [400 450 500 550 600 650 700 725]; % wavelenght


%for all wavelength read file and calcuate PD and DO
for i = 1:size(lambda,2)
    mcml_data_d = Readmcml("data_files/outputs/sample_d_n_" + lambda(i) + ".mco");
    mcml_data_s = Readmcml("data_files/outputs/sample_s_n_" + lambda(i) + ".mco");
    [PD(i),DO(i)] = pd_do(mcml_data_d,mcml_data_s,PLOTON);
end
%% 

%make figure with all skin Layers and the PD DO plot
figure
hold on
Fz_size = size(mcml_data_d.Fz,1);
PD = PD.*mcml_data_d.dz;
plot(lambda,PD);
DO = DO.*mcml_data_d.dz;
plot(lambda,DO);
depth = 0;
for i = 1:size(mcml_data_d.d)
    depth = depth + mcml_data_d.d(i); 
    plot([lambda(1) lambda(end)],[depth depth]);
end
xlabel('wavelength [nm]')
ylabel('PD')
axis([lambda(1) lambda(end) 0 Fz_size*mcml_data_d.dz])
set(gca, 'YDir','reverse');
%%

function [PD DO] = pd_do(mcml_data_d,mcml_data_s,PLOTON)
%%%%
% calculation for PD and DO
% PD = 63.2% of the area under the diastolic or systolic flow curves
% DO = 63.2% of the area under the difference of diastolic-systolic flow curves.
%%%%

    pulsation_pattern_w_comp = [0 0 0 2/3 18/3 1/3];
    pulsation_pattern_w_ref = [0 0 1/3 2/3 1 1/3];

    %determine PD
    Fz = mcml_data_d.Fz;
    PD_threshold = 0.632* sum(Fz);
    F_pd = zeros(size(Fz));
    F_pd(1) = sum(Fz);
    first = 1;
    for depth = 2:size(Fz)
        F_pd(depth) = F_pd(depth-1)-Fz(depth);
        if F_pd(depth) <= PD_threshold && first == 1
            PD = depth;
            first = 0;
        end
    end

    %determine the sum function for Fz systolic
    Fz_s = mcml_data_s.Fz;
    F_pd_s = zeros(size(Fz_s));
    F_pd_s(1) = sum(Fz_s);
    for depth = 2:size(Fz_s)
        F_pd_s(depth) = F_pd_s(depth-1)-Fz_s(depth);
    end

    %determine Delta Fz
    delta_Fz = zeros(size(Fz_s));
    layer = 1;
    sum_layer = mcml_data_d.d(1);
    for depth = 1:size(Fz)
        delta_Fz(depth) = abs((F_pd(depth) - F_pd_s(depth)) * pulsation_pattern_w_ref(layer));
        if depth*mcml_data_d.dz > sum_layer
            layer = layer+1;
            sum_layer = sum_layer+mcml_data_d.d(layer);
        end
    end

    %determine PD
    DO = 0;
    DO_threshold = 0.632* sum(delta_Fz);
    delta_sum = 0;
    for depth = 1:size(delta_Fz)
        delta_sum = delta_sum + delta_Fz(depth);
        if delta_sum >= DO_threshold
            DO = depth;
            break;
        end
    end


    %Plot
    if PLOTON
        f   = linspace(0,1,size(Fz,1));
        f_delta = linspace(0,1,size(delta_Fz,1));
        F_pd = F_pd.';
        figure
        F_pd_plot = F_pd./sum(Fz);
        PD_plot = PD/size(Fz,1);
        DO_plot = DO/size(delta_Fz,1);
        DO_max = max(delta_Fz);
        delta_Fz = delta_Fz./DO_max;
        plot(F_pd_plot,f,[0 1],[PD_plot PD_plot],[F_pd_plot(PD) F_pd_plot(PD)],[0 1],delta_Fz,f_delta,[0 1],[DO_plot DO_plot]);
        xlabel('Flux in percent')
        ylabel('PD in cm')
        set(gca, 'YDir','reverse');
        grid;
    end
end
