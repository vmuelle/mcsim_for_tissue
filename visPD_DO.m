close all;
clear all;



%%%%
%This file visualizes the penetretion-depth (PD) and depth-origin (DO) of
%light in skin tissue. The outputfiles of MCML are used for this. 
%%%%

PLOTON = 0;%If PLOTON = 1, the Flux curves are shown for each input file. If PLOTON = 0 only the final plot over all Wavelength is shown.
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
A_size = size(mcml_data_d.Az,1);
A_size*mcml_data_d.dz;
plot(lambda,PD);
plot(lambda,DO);
depth = 0;
for i = 1:size(mcml_data_d.d)
    depth = depth + mcml_data_d.d(i); 
    plot([lambda(1) lambda(end)],[depth depth]);
end
xlabel('wavelength [nm]')
ylabel('PD')
axis([lambda(1) lambda(end) 0 A_size*mcml_data_d.dz])
set(gca, 'YDir','reverse');
%%

function [PD DO] = pd_do(mcml_data_d,mcml_data_s,PLOTON)
%%%%
% calculation for PD and DO
% PD = 63.2% of the area under the diastolic or systolic flow curves
% DO = 63.2% of the area under the difference of diastolic-systolic flow curves.
%%%%
    Az = mcml_data_d.Az;
    A_sum = sum(Az,1);
    A_size = size(Az,1);
    A_flux = zeros(1,A_size);
    A_flux(1) = A_sum;
    PD_sum = Az(1);
    for depth = 2:A_size
        A_flux(depth) = A_flux(depth-1)-Az(depth);
        if(PD_sum > 0.632*A_sum)
            PD = depth/A_size;
            break;
        end
        PD_sum = PD_sum + Az(depth);
    end

    delta_d_s_Az = abs(mcml_data_d.Az - mcml_data_s.Az);
    delta_d_s_sum = sum(delta_d_s_Az);
    DO_sum = 0;
    for depth = 1:A_size
        if(DO_sum > 0.632*delta_d_s_sum)
            DO = depth/A_size;
            break;
        end
        DO_sum = DO_sum + delta_d_s_Az(depth);
    end

    if PLOTON
        A_flux = A_flux./A_sum; 
        f   = [1:-1/(A_size-1):0];
        figure
        plot(f,1-A_flux,[0 1],[PD PD],f,delta_d_s_Az,[0 1],[DO DO]);
        %plot(f,1-A_flux,[0 1],[PD PD]);
        xlabel('Flux in percent')
        ylabel('PD in cm')
        set(gca, 'YDir','reverse');
        grid;
    end
end
