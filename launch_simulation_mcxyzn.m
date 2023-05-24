function launch_simulation_mcxyzn(varargin)
    if(nargin == 0) 
        p = 0.3;
        make_circle_ = 0;
        radius = 0;
        center = 0;
        tissue_type.mua = [0,0,0,0,0,0,0,0];
    elseif (nargin == 1)
        p = varargin{1};
        make_circle_ = 0;
        radius = 0;
        center = 0;
        tissue_type.mua = [0,0,0,0,0,0,0,0];
    elseif  (nargin == 4)
        p = varargin{1};
        %radius = 50/2;
        %center = [200,200,800]/2;
        %tissue_type = 3;
        radius = varargin{2};
        center = varargin{3};
        tissue_type = varargin{4};
        make_circle_ = 1;
    end
    
    
    %Transform parameters into sufficient values
    tissue = Tissue();
    coeff = tissue.makeTissueList('diastolic','normal',p);
    %coeff = coeff(1,:);
    for i = 1:size(coeff,1)
        simulation_per_freq(coeff(i,:),'d','n',radius,center,tissue_type.mua(i),make_circle_);
    end
    coeff = tissue.makeTissueList('diastolic','compressed',p);
    %coeff = coeff(1,:);
    for i = 1:size(coeff,1)
        simulation_per_freq(coeff(i,:),'d','c',radius,center,tissue_type.mua(i),make_circle_);
    end
    coeff = tissue.makeTissueList('systolic','normal',p);
    %coeff = coeff(1,:);
    for i = 1:size(coeff,1)
        simulation_per_freq(coeff(i,:),'s','n',radius,center,tissue_type.mua(i),make_circle_);
    end
    coeff = tissue.makeTissueList('systolic','compressed',p);
    %coeff = coeff(1,:);
    for i = 1:size(coeff,1)
        simulation_per_freq(coeff(i,:),'s','c',radius,center,tissue_type.mua(i),make_circle_);
    end
end
    
function simulation_per_freq(coeff,d_s,c_n,radius,center,tissue_type,make_circle_)

    %Parameters for the simulation
    %photons = 100000;
    dz = 0.001; 
    n_dr = 600; %1000;
    %n_dz = 800; 
    n_dz = 0;
    for i = 1:size(coeff,2)
        n_dz = n_dz + coeff(i).d/dz;
    end
    
    %%% Basic MC configuration %%%
    cfg.SAVEON = 1; % 1 = save myname_T.bin, myname_H.mci 
                    % 0 = don't save. Just check the program.
    cfg.name = "data_files\outputs\mcxyzn\moco_params_"+d_s+"_"+c_n+"_"+coeff(1).nm;
    %cfg.name = "moco_params_d_n_"+coeff(1).nm;
    cfg.time = 1;               %Simulation time in min
    
    cfg.binsize = dz;        %Length of a voxel
   
    n_dz = cast(n_dz,'uint32');
    n_dz = cast(n_dz,'double');
    cfg.dim = [n_dr,n_dr,n_dz]; %Number of voxels in each direction [Nx,Ny,Nz]
    
    
    cfg.mcflag       = 3;   % launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt. 
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
    cfg.radius   = 0.0001;         % 1/e radius of beam at tissue surface
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
        cfg.gv   = [cfg.gv coeff(layer).g];
        cfg.nv   = [cfg.nv coeff(layer).n];
    end
    
    %%% Setting up simulation volume %%% 
    T = double(ones(cfg.dim(1),cfg.dim(2),cfg.dim(3))); 
    zsurf = 0.0100;  % position of air/skin surface
    
    
    cfg.muav = [cfg.muav tissue_type];
    cfg.musv = [cfg.musv coeff(1).mus];
    cfg.gv   = [cfg.gv coeff(1).g];
    cfg.nv   = [cfg.nv coeff(1).n];

    T_t = make_tissue(cfg,T,coeff,dz);

    if (make_circle_)
        T_t = make_circle(cfg,T_t,radius,center,size(cfg.muav,2));
        cfg.name = "data_files\outputs\mcxyzn\moco_params_circle_"+d_s+"_"+c_n+"_"+coeff(1).nm;
    end

    cfg.T = T_t;
    
    disp('Tissue finished')

    clear T
    clear T_t;
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
                        if(lay ~= length(coeff))    
                            lay = lay+1;
                        end
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
                 %if ((round(sqrt(double((center(1)-i)^2+(center(2)-j)^2+(center(3)-iz)^2)))<radius_c))
                 if ((round(sqrt(double((center(2)-j)^2+(center(3)-iz)^2)))<radius_c))
                    T(i,j,iz) = tissue_type;
                 end
             end
         end
     end % iz
end
 



