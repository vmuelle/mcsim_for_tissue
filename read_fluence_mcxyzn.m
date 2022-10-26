clear all;
close all;

freq = [470,530,660,770,810,940,1020,1050];
j = 1;
PD = zeros(length(freq));
for i = freq
    Fz = getfluence(i);
    PD(j) = getPD(Fz);
    j=j+1;
end
make_plot(PD,freq);


function PD = getPD(Fz)
    PD_threshold = 0.632* sum(Fz);
    Fz_sum = 0;
    for depth = 1:size(Fz,2)
        Fz_sum = Fz_sum+Fz(depth);
        if Fz_sum >= PD_threshold
            PD = depth;
            break;
        end
    end
end

function make_plot(PD,lambda)
    figure
    mcml_data_d = Readmcml("data_files/outputs/sample_d_n_470.mco");
    depth = 0;
    for i = 1:size(mcml_data_d.d)
        depth = depth + mcml_data_d.d(i); 
        name = sprintf('Layer %d', i);
        plot([lambda(1) lambda(end)],[depth depth],'DisplayName',name);
        hold on
    end
    PD = PD./800;
    plot(lambda,PD)
    set(gca, 'YDir','reverse');
end



function Fz = getfluence(i)
    name = strcat('data_files/outputs/mcxyzn/moco_params_d_n_',num2str(i,'%d'));
    
    filename = sprintf('%s_H.mci',name);
    disp(['loading ' filename])
    fid = fopen(filename, 'r');
    A = fscanf(fid,'%f',[1 Inf])';
    fclose(fid);
    
    Nx = A(2);
    Ny = A(3);
    Nz = A(4);
    
    %% Load Fluence rate F(x,y,z) 
    filename = sprintf('%s_F.bin',name);
    disp(['loading ' filename])
    tic
        fid = fopen(filename, 'rb');
        [Data count] = fread(fid, Nx*Ny*Nz, 'float');
        fclose(fid);
    toc
    F = reshape(Data,Nx,Ny,Nz);
    
    for i = 1:size(F,3)
        Fz(i) = sum(sum(F(:,:,i)));
    end
    x = [0:1:800-1];
    figure
    plot(x,Fz);
end