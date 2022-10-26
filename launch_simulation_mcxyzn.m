function launch_simulation_mcxyzn()
    clear all;
    
    %Transform parameters into sufficient values
    tissue = Tissue();
    coeff = tissue.makeTissueList('diastolic','normal');
    %coeff = coeff(1,:);
    for i = 1:size(coeff,1)
        simulation_per_freq(coeff(i,:));
    end
end
    
function simulation_per_freq(coeff)

    %Parameters for the simulation
    photons = 100000;
    dz = 0.001; 
    n_dr = 200; %1000;
    n_dz = 800;
    
    %%% Basic MC configuration %%%
    cfg.SAVEON = 1; % 1 = save myname_T.bin, myname_H.mci 
                    % 0 = don't save. Just check the program.
    cfg.name = "data_files\outputs\mcxyzn\moco_params_d_n_"+coeff(1).nm;
    %cfg.name = "moco_params_d_n_"+coeff(1).nm;
    cfg.time = 1;               %Simulation time in min
    
    cfg.binsize = dz;        %Length of a voxel
    sum_d = 0;
    for i = 1:size(coeff,2)
        sum_d = sum_d + coeff(i).d;
    end
    Nz = sum_d/dz;
    Nz = cast(Nz,'uint32');
    cfg.dim = [n_dr,n_dr,n_dz]; %Number of voxels in each direction [Nx,Ny,Nz]
    
    
    cfg.mcflag       = 0;   % launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt. 
                            % 3 = rectangular beam (use xfocus,yfocus for x,y halfwidths)
    cfg.launchflag   = 0;   % 0 = let mcxyz.c calculate launch trajectory
                            % 1 = manually set launch vector.
    cfg.boundaryflag = 1;   % 0 = no boundaries, 1 = escape at boundaries
                            % 2 = escape at surface only. No x,y, bottom z
                            % boundaries
    cfg.gradientflag = 1;   % 0 = Fresnel's law using the voxel faces [inaccurate for curved and oblique geometries]
                            % 1 = gradient-based Fresnel's laws w. smoothing
                            % 2 = gradient-based Fresnel's laws w. smoothing and interpolation 
    cfg.srcpos   = [0,0,0.0001];   % Set position of source [x,y,z];
    cfg.srcfocus = [0,0,inf]; % Set position of focus, so mcxyz can calculate launch trajectory
                              % [xfocus,yfocus,zfocus];
    cfg.radius   = 0.15;         % 1/e radius of beam at tissue surface
    cfg.waist    = 0.15;         % 1/e radius of beam at focus
    
    % only used if launchflag == 1 (manually set launch trajectory) / ux^2 + uy^2 + uz^2 = 1
    cfg.launchvec = [0.7,0.4,sqrt(1 - 0.7^2 - 0.4^2)]; 
    
    cfg.Nt   = size(coeff,2);
    cfg.muav = [];
    cfg.musv = [];
    cfg.gv = [];
    cfg.nv = [];
    for layer = 1:(cfg.Nt)    
        cfg.muav = [cfg.muav coeff(layer).mua];
        cfg.musv = [cfg.musv coeff(layer).mus];
        cfg.gv   = [cfg.muav coeff(layer).g];
        cfg.nv   = [cfg.muav coeff(layer).n];
    end
    
    %%% Setting up simulation volume %%% 
    T = double(ones(cfg.dim(1),cfg.dim(2),cfg.dim(3))); 
    zsurf = 0.0100;  % position of air/skin surface
    
    
    T_t = make_tissue(cfg,T,coeff,dz);
    
    radius = 50/2;
    center = [200,200,800]/2;
    %T_tc = make_circle(cfg,T_t,radius,center,3);
    
    cfg.T = T_t;
    
    disp('Tissue finished')
    create_simfiles(cfg);
    %%
    launch_simulation(cfg);
    %[fluence] = plot_mcresults(cfg);

end

function T= make_tissue(cfg,T_,coeff,dz)
    T = T_;
    lay = 1;
    for iz=1:cfg.dim(3) % for every depth z(iz)
        for i = 1:cfg.dim(2)
            for j = 1:cfg.dim(1)
                if(lay < cfg.Nt+1)
                    sum_d = 0;
                    for k = 1:lay
                        sum_d = sum_d + coeff(k).d;
                    end
                    if(iz<cast(sum_d/dz,'uint32'))
                        T(i,j,iz) = lay;
                    else 
                        lay = lay+1;
                        T(i,j,iz) = lay;
                    end
                end
            end
        end
    end
end


%einfügen eines zusätzlichen blutgefäß
function T = make_circle(cfg,T_,radius_c,center,tissue_type)
    T = T_;
    for iz=1:cfg.dim(3)
         for i = 1:cfg.dim(2)
             for j = 1:cfg.dim(1)
                 if ((round(sqrt(double((center(1)-i)^2+(center(2)-j)^2+(center(3)-iz)^2)))<radius_c))
                    T(i,j,iz) = tissue_type;
                 end
             end
         end
     end % iz
end
 



