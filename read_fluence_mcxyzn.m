function read_fluence_mcxyzn(PLOTON)
    freq = [470,530,660,770,810,940,1020,1050];
    j = 1;
    PD = zeros(length(freq),1);
    for i = freq
        Fz = getfluence(i,'n',PLOTON);
        PD(j) = getPD(Fz);
        j=j+1;
    end
    make_plot(PD,freq,'n',length(Fz));

    j = 1;
    PD = zeros(length(freq),1);
    for i = freq
        Fz = getfluence(i,'c',PLOTON);
        PD(j) = getPD(Fz);
        j=j+1;
    end
    make_plot(PD,freq,'c',length(Fz));

    
end

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

function make_plot(PD,lambda,n_c,len)
    figure
    hold on
    filename = "data_files/outputs/sample_d_"+n_c+"_470.mco";
    mcml_data_d = Readmcml(filename);
    Fz_size = size(mcml_data_d.Fz,1);
    PD = PD.*mcml_data_d.dz;
    plot(lambda,PD,'DisplayName','PD')

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
    if(length(mcml_data_d.d) == 6)
        title("Penetration Depth for normal skin mcxyzn")
    else
        title("Penetration Depth for compressed skin mcxyzn")
    end
    hold off
    legend

end



function Fz = getfluence(i,n_c,PLOTON)
    name = strcat('data_files/outputs/mcxyzn/moco_params_d_',n_c,'_',num2str(i,'%d'));
    
    filename = sprintf('%s_H.mci',name);
    %disp(['loading ' filename])
    fid = fopen(filename, 'r');
    A = fscanf(fid,'%f',[1 Inf])';
    fclose(fid);
    
    Nx = A(2);
    Ny = A(3);
    Nz = A(4);
    
    %% Load Fluence rate F(x,y,z) 
    filename = sprintf('%s_F.bin',name);
    %disp(['loading ' filename])
    tic
        fid = fopen(filename, 'rb');
        [Data count] = fread(fid, Nx*Ny*Nz, 'float');
        fclose(fid);
    toc
    F = reshape(Data,Nx,Ny,Nz);
    
    for i = 1:size(F,3)
        Fz(i) = sum(sum(F(:,:,i)));
    end
    x = [0:1:length(Fz)-1];
    %x = [0:1:800-1];
    if(PLOTON)
    figure
    plot(x,Fz);
    end
end