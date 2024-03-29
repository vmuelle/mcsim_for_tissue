%close all;
%clear all;

function read_fluence_mcml_chatterjee(PLOTON)

%%%%
%This file visualizes the penetretion-depth (PD) and depth-origin (DO) of
%light in skin tissue. The outputfiles of MCML are used for this. 
%%%%

%PLOTON = 1;%If PLOTON = 1, the Flux curves are shown for each input file. If PLOTON = 0 only the final plot over all Wavelength is shown.
PD = zeros(1,8); % Parameter for penetration-depth
%lambda = [400 450 500 550 600 650 700 750 800 850 900 950 1000]; % wavelenght
lambda = [470 530 660 770 810 940 1020 1050];
%lambda = [400 450 500 550 600 650 700 725];
PD_c = zeros(1,8); % Parameter for penetration-depth
DO_c = zeros(1,8); % Parameter for depth-origin 

%for all wavelength read file and calcuate PD and DO
for i = 1:size(lambda,2)
    mcml_data = Readmcml("data_files/outputs/chatterjee_" + lambda(i) + ".mco");
    PD(i) = pd_do(mcml_data,PLOTON,lambda(i));
end
make_figure(PD,mcml_data,lambda);

end


function make_figure(PD,mcml_data_d,lambda)
%make figure with all skin Layers and the PD DO plot
figure
hold on
Fz_size = size(mcml_data_d.Fz,1);
PD = PD.*mcml_data_d.dz;
plot(lambda,PD,'DisplayName','PD');
depth = 0;
for i = 1:size(mcml_data_d.d)
    depth = depth + mcml_data_d.d(i); 
    name = sprintf('Layer %d', i);
    plot([lambda(1) lambda(end)],[depth depth],'DisplayName',name);
end
xlabel('wavelength [nm]')
ylabel('skin depth [mm]')
axis([lambda(1) lambda(end) 0 Fz_size*mcml_data_d.dz])
set(gca, 'YDir','reverse');
title("Penetration Depth for normal skin")
hold off
legend
end


%%

function PD  = pd_do(mcml_data_d,PLOTON,lambda)
%%%%
% calculation for PD and DO
% PD = 63.2% of the area under the diastolic or systolic flow curves
% DO = 63.2% of the area under the difference of diastolic-systolic flow curves.
%%%%

    pulsation_pattern_w_comp = [0 0 0 2/3 18/3 1/3];
    pulsation_pattern_w_ref = [0 0 1/3 2/3 1 1/3];

    %determine PD
    %Fz = mcml_data_d.Fz;
    %PD_threshold = 0.632* sum(Fz);
    %F_pd = zeros(size(Fz));
    %F_pd(1) = sum(Fz);
    %first = 1;
    %for depth = 2:size(Fz)
    %    F_pd(depth) = F_pd(depth-1)-Fz(depth);
    %    if F_pd(depth) <= PD_threshold && first == 1
    %        PD = depth;
    %        first = 0;
    %    end
    %end
    Fz = mcml_data_d.Fz;
    layer_boundarys = mcml_data_d.d./mcml_data_d.dz;
    for i = 2:size(layer_boundarys)
        layer_boundarys(i) = layer_boundarys(i) +layer_boundarys(i-1);
    end
    for i = 1:length(layer_boundarys)-1
        %Fz(layer_boundarys(i))
        %Fz(layer_boundarys(i)-1)
        if Fz(layer_boundarys(i)) > Fz(layer_boundarys(i)-1)
            difference = max(Fz(layer_boundarys(i):layer_boundarys(i)+10))-Fz(layer_boundarys(i)-1);
            Fz(layer_boundarys(i):end) = Fz(layer_boundarys(i):end)-difference;
        end
    end




    %Fz = Fz-min(Fz);
    %Fz = Fz./max(Fz);
    PD_threshold = 0.632* sum(Fz);
    Fz_sum = 0;
    for depth = 1:size(Fz)
        Fz_sum = Fz_sum+Fz(depth);
        if Fz_sum >= PD_threshold
            PD = depth;
            break;
        end
    end


    %figure
    %pp   = linspace(0,1,size(Fz,1));
    %plot(pp,Fz.')
    %figure
    %plot(pp,F_pd)

    %determine the sum function for Fz systolic
    %Fz_s = mcml_data_s.Fz;
    %F_pd_s = zeros(size(Fz_s));
    %F_pd_s(1) = sum(Fz_s);
    %for depth = 2:size(Fz_s)
    %    F_pd_s(depth) = F_pd_s(depth-1)-Fz_s(depth);
    %end




    %Plot
    if PLOTON
        f   = linspace(0,size(Fz,1),size(Fz,1));
        figure
        F_pd_plot = Fz;%/max(Fz);
        %plot(f,F_pd_plot,[PD PD],[0 1])
        %PD_plot = PD/size(Fz,1);
        %DO_plot = DO/size(delta_Fz,1);
        hold on
        plot(f,F_pd_plot,'DisplayName','Flux in percent')
        plot([PD PD],[0 1],'DisplayName','PD(at 63 percent)')
        %plot([F_pd_plot(PD) F_pd_plot(PD)],[0 1],'DisplayName','')
        %plot(f,delta_Fz,'DisplayName','Delta Flux')
        %plot([DO DO],[0 1],'DisplayName','DO');
        hold off 
        legend
        xlabel('depth (grid elements)')
        ylabel('Flux F/max(F)')
        %set(gca, 'YDir','reverse');
        %view([90 -90])
        title(sprintf("Flux and Delta flux for normal skin at %d nm",lambda))
        grid;
    end
end
