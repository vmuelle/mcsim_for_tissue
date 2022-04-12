%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculation for the coefficients of absorption and scattering for 
% different skin tissues to input into mcml and execution by mcml
%
% Autor:         Viktor Müller
% Datum:         09.04.22
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global c_w_0 spO2 svO2 p anisotrophy;
global absorption_hbO2 absorption_hb absorption_water absorption_fat c_b c_w c_f v_d sc ratio_a_nv nm d;

%get the parameters for tissue with different properties
tissue_list_dn = makeTissueList('diastolic','normal');
tissue_list_dc = makeTissueList('diastolic','compressed');
tissue_list_sn = makeTissueList('systolic','normal');
tissue_list_sc = makeTissueList('systolic','compressed');


%set parameters for simulation 
photons = 100000;             %no of photons
dz = 0.001;                   %separation between grid lines (cm) in z direction                                                                                                                                                                 
dr = 0.001;                   %separation between grid lines (cm) in r direction
n_dr = 1000;                  %no of grid elements in r direction        
n_da = 90;                    %no of parts of an 90 degree angle, that is spanned between the photon exiting direction and the surface normal
sim = [photons,dz,dr,n_dr,n_da];

%print mcml input files and launch simulation
filename = 'data_files/inputs/coefficients_d_n.mci';
fileID = fopen(filename,'w');
print_to_file(tissue_list_dn,fileID,'d','n',sim);
fclose(fileID);
launch_simulation(filename);

filename = 'data_files/inputs/coefficients_d_c.mci';
fileID = fopen(filename,'w');
print_to_file(tissue_list_dc,fileID,'d','c',sim);
fclose(fileID);
launch_simulation(filename);

filename = 'data_files/inputs/coefficients_s_n.mci';
fileID = fopen(filename,'w');
print_to_file(tissue_list_sn,fileID,'s','n',sim);
fclose(fileID);
launch_simulation(filename);

filename = 'data_files/inputs/coefficients_s_c.mci';
fileID = fopen(filename,'w');
print_to_file(tissue_list_sc,fileID,'s','c',sim);
fclose(fileID);
launch_simulation(filename);

function launch_simulation(filename)
%function launch_simulation(filename)
% launches the simulation via mcml with the file given by 'filename'
    if (ispc)
        bin_name = 'mcml.exe';
    else
        fprintf('Could not find the appropriate binary \n');
    end
    filepath = which(bin_name);
    system_command_string = [sprintf('%s',filepath),' ',filename];

    [status] = system(system_command_string);
end

function tissueList = makeTissueList(d_,n_)
%function tissue = makeTissueList(d_,n_)
%  Returns a list of tissue structures with wavelengths from 400 to 1000 nm
%  with 50 nm steps (so it's 13 wavelenths) with properties d_: 'diastolic' or 'systolic and n_: 'normal' or
%  'compressed'
    j = 1;
    for i = 400:50:700  %for i = 400:50:1000 
        tissue = makeTissue(i,d_,n_);
        tissue(1).nm = i;
        tissueList(j,:) = tissue;
        j=j+1;
    end
end

function tissue = makeTissue(nm_,d_,n_)
%function tissue = makeTissue(nm_,d_,n_)
%   Returns the tissue optical properties at the wavelength nm with
%   properties d_: 'diastolic' or 'systolic and n_: 'normal' or
%   'compressed'
%       tissue = [name; n; mua; mus; g; d]';
%       uses SpectralLIB.mat
    global c_w_0 spO2 svO2 p anisotrophy sc_b;
    global absorption_hbO2 absorption_hb absorption_water absorption_fat c_b c_w c_f v_d sc ratio_a_nv nm d;
      
    nm = nm_;           %wavelenth
    c_w_0 = 0.65;       %coefficient that accounts for background measurement
    spO2 = 0.97;        %aterial oxygen saturation
    svO2 = 0.67;        %venious oxgen saturation
    p = 0.1;            %pulsation increase
    anisotrophy = 0.9;  %anisotrophy
    sc_b = 0.1;         %decaying factor for scattering
    if(nm>580)                                                                          
        sc_b = 0.05;
    end

    % Load spectral library
    load spectralLIBv2.mat %load spectralLIB.mat
    %   muadeoxy      701x1              5608  double              
    %   muamel        701x1              5608  double              
    %   muaoxy        701x1              5608  double              
    %   muawater      701x1              5608  double              
    %   musp          701x1              5608  double              
    %   nmLIB         701x1              5608  double              
    absorption_hbO2 = interp1(nmLIB,muaoxy,nm);     %absorption coefficient hb02
    absorption_hb = interp1(nmLIB,muadeoxy,nm);     %absorption coefficient hb
    absorption_water = interp1(nmLIB,muawater,nm);  %absorption coefficient water
    absorption_fat = interp1(nmLIB,muafat,nm);      %absorption coefficient fat
    
    
    

    if(strcmp(n_,'normal'))
        %fixed parameters given by moco per layer
        n = [1.33 1.37 1.40 1.40 1.40 1.44];    %reflective index
        d = [0.08 0.015 0.008 0.12 0.05 0.5];   %thickness  
        c_b = [0 0.004 0.02 0.004 0.04 0.03];   %concentration of blood
        c_w = [0.2 0.65 0.65 0.65 0.65 0.05];   %concentration of water
        c_f = [0 0 0 0 0 0.4];                  %concentration of fat
        v_d = [0 10 20 20 40 50];               %vessel diameters
        sc = [15 20 20 20 20 10];               %scattering calibration constant
        ratio_a_nv = [0 0.5 0.5 0.5 0.5 0.5];   %ratio aterial blood (not venious)

        tissue(1).name = 'EPIDERMIS';
        tissue(2).name = 'CAPILLARY LOOPS';
        tissue(3).name = 'UPPER PLEXUS';
        tissue(4).name = 'RETICULAR PLEXUS';
        tissue(5).name = 'DEEP PLEXUS';
        tissue(6).name = 'HYPODERMIS';
        
        for i = 1:6
            tissue(i).n = n(i);
            tissue(i).d = d(i);
            tissue(i).g = anisotrophy;
            tissue(i).mus = scattering(i);
        end
        
         tissue(1).mua   = absorption_first_layer();

        for i = 2:6
            if(strcmp(d_,'diastolic'))
                tissue(i).mua = absorption_diastolic_state(i);
            elseif(strcmp(d_,'systolic'))
                tissue(i).mua = absorption_systolic_state(i);
            end
        end

    elseif(strcmp(n_,'compressed'))
        %fixed paraters, given by moco for compressed skin
        n = [1.33 1.37 1.40 1.40 1.44];         %reflective index
        d = [0.08 0.008 0.004 0.1 0.2];         %thickness
        c_b = [0 0.0012 0.0024 0.024 0.036];    %concentration of blood
        c_w = [0.05 0.15 0.15 0.15 0.35];       %concentration of water
        c_f = [0 0 0 0 0.4];                    %concentration of fat
        v_d = [0 10 20 20 40];                  %vessel diameters
        sc = [15 20 20 20 10];                  %scattering calibration constant
        ratio_a_nv = [0 1 1 1 0.75];            %ratio aterial blood (not venious)

        tissue(1).name = 'EPIDERMIS';
        tissue(2).name = 'CAPILLARY LOOPS';
        tissue(3).name = 'UPPER PLEXUS';
        tissue(4).name = 'RETICULAR PLEXUS & DEEP PLEXUS';
        tissue(5).name = 'HYPODERMIS';

        for i = 1:5
            tissue(i).n = n(i);
            tissue(i).d = d(i);
            tissue(i).g = anisotrophy;
            tissue(i).mus = scattering(i);
        end
        
        tissue(1).mua   = absorption_first_layer();

        for i = 2:5
            if(strcmp(d_,'diastolic'))
                tissue(i).mua = absorption_diastolic_state(i);
            elseif(strcmp(d_,'systolic'))
                tissue(i).mua = absorption_systolic_state(i);
            end
        end
    end
end

function [f_a, f_v, c_b_corrected] = blood_concentration_effective_diastolic(layer)
%function [f_a, f_v, c_b_corrected] = blood_concentration_effective_diastolic(layer)
% Returns the fraction of aterial/venious blood 'f_a'/'f_v' and the effective
% blood concentration 'c_b_corrected' for a given layer in diastolc state
    global spO2 svO2 absorption_hbO2 absorption_hb c_b v_d ratio_a_nv;

    f_a = ratio_a_nv(layer) * c_b(layer) * arterial_venous_c_b(((1-spO2)*absorption_hb)+(spO2*absorption_hbO2),v_d(layer));
    f_v = (1-ratio_a_nv(layer)) * c_b(layer) * arterial_venous_c_b(((1-svO2)*absorption_hb)+(svO2*absorption_hbO2),v_d(layer));
    c_b_corrected = f_a + f_v;
end


function [d_s, f_a_s, f_v, c_b_corrected] = blood_concentration_effective_systolic(layer)
%function [d_s, f_a_s, f_v, c_b_corrected] = blood_concentration_effective_systolic(layer)
% Returns the corrected thickness 'd_s', the fraction of aterial/venious blood 'f_a_s'/'f_v' and the effective
% blood concentration 'c_b_corrected' for a given layer in systolc state
    global spO2 p absorption_hbO2 absorption_hb v_d d;

    [f_a, f_v, c_b_corrected] = blood_concentration_effective_diastolic(layer);

    f_a_s = f_a + p*arterial_venous_c_b(((1-spO2)*absorption_hb)+(spO2*absorption_hbO2),v_d(layer));
    c_b_corrected = c_b_corrected+p*arterial_venous_c_b(((1-spO2)*absorption_hb)+(spO2*absorption_hbO2),v_d(layer));
    d_s = d(layer)*(1+c_b_corrected*p);
end



function afl =  absorption_first_layer()
%function afl =  absorption_first_layer()
% Returns the absorption coefficient 'afl' for first layer (epidermis)
    global c_b c_w absorption_water nm;

    afl = (c_w(1)*absorption_water) + (1-c_b(1)) * 0.5*(0.244+(85.3*exp((-nm-154)/66.2)));
end
        
       
function a_base = absorption_base(layer)
%function a_base = absorption_base(layer)
% Returns the base absorption 'a_base' for a given 'layer'
    global c_w_0 c_w nm;
   
    factor = 0.5;
    if(layer == length(c_w))
        factor = 0.25;
    end
    a_base = factor * (c_w(layer)/c_w_0)* (0.244 + 16.82* exp((-nm-400)/80.5));
end


function a_c = absorption_without_blood(layer)
%function a_c = absorption_without_blood(layer)
% Returns the absorption coefficient 'a_c' for a given 'layer' WITHOUT blood
    global absorption_fat c_w c_f absorption_water;

    a_c = c_f(layer)*absorption_fat + (1-c_f(layer))* c_w(layer) * absorption_water + (1-c_f(layer))*(1-c_w(layer))*absorption_base(layer);
end


function mua = absorption_diastolic_state(layer)
%function mua = absorption_diastolic_state(layer)
% Returns the absorption coefficient 'mua' for a given 'layer' WITH blood in diastolic state
    global spO2 svO2 absorption_hbO2 absorption_hb c_b absorption_water

    [f_a, f_v, c_b_corrected] = blood_concentration_effective_diastolic(layer);
    mua = f_a*((1-spO2)*absorption_hb+spO2*absorption_hbO2) ...
      + f_v*((1-svO2)*absorption_hb+svO2*absorption_hbO2) ...
      + c_b_corrected*absorption_water+(1-c_b(layer))*absorption_without_blood(layer);
end


function mua = absorption_systolic_state(layer)
%function mua = absorption_systolic_state(layer)
% Returns the absorption coefficient 'mua' for a given 'layer' WITH blood in systolic state
    global spO2 svO2 absorption_hbO2 absorption_hb absorption_water d;

    [d_s, f_a_s, f_v, c_b_corrected] = blood_concentration_effective_systolic(layer);
    mua = d(layer)/d_s ...
        * (... 
        f_a_s*((1-spO2)*absorption_hb+spO2*absorption_hbO2) ...
        + f_v*((1-svO2)*absorption_hb+svO2*absorption_hbO2)...
        + c_b_corrected*absorption_water+(1-c_b_corrected)*absorption_without_blood(layer)...
        );
end


function mus = scattering(layer)
%function mus = scattering(layer)
% Returns the scattering coefficient 'mus' for a given 'layer'
    global sc sc_b nm;

    mus = sc(layer)* nm^(-sc_b);
end

function avcb = arterial_venous_c_b(absorption,v_d)
%function avcb = arterial_venous_c_b(absorption,v_d)
% Returns a correction factor 'avcb' that compensates the vessel walls to
% get an apperent blood volume.
% property 'absorption' hast to be a blood absorption parameter. 
% property 'v_d' is a vessel diameter
    avcb = 1/(1+1.007*(absorption*v_d/2)^1.228);
end



function print_to_file(tissueList,fileID,d_,n_,sim)
%function print_to_file(tissueList,fileID,d_,n_,sim)
% prints the properties of the 'tissueList' into a .mci formatted file with name 'fileID', that can
% be used by the mcml simulation software.
% sim = [photons,dz,dr,n_dr,n_da];
    number_of_runs = size(tissueList,1);
    fprintf(fileID,['' ...
        '##########################################\n'...
        '# coefficients_%c_%c.mci\n'...
        '# 	This is an input files for MCM\n'...
        '#	Lengths are in cm, mua and mus are in 1/cm.\n'...
        '#\n'...
        '#	Multiple runs are stipulated.\n'...
        '##########################################\n'...
        '\n'...
        '1.0                      	# file version\n'...
        '%d                        	# number of runs\n'],d_,n_,number_of_runs);
    for freq = 1:number_of_runs
        tissue = tissueList(freq,:);
        n_dz = 0;
        for layer = 1:size(tissue,2)
            n_dz = n_dz + tissue(layer).d;
        end
        n_dz = n_dz/sim(2);
        n_dz = cast(n_dz,'uint32');
        fprintf(fileID,[''...
            '#### SPECIFY DATA FOR RUN %d at lambda = %d\n'...
            '#InParm                    	# Input parameters. cm is used.\n'...
            'data_files/outputs/sample_%c_%c_%d.mco 	A            	# output file name, ASCII.\n'...
            '%d	                  	    # No. of photons\n'...
            '%f	%f               	# dz, dr [cm]\n'...
            '%d	%d	%d                	# No. of dz, dr, da.\n' ...
            '\n'...
            '%d                        	# Number of layers\n'...
            '#n	mua	mus	g	d         	# One line for each layer\n'...
            '1                         	# n for medium above\n'], ...
            freq,tissue(1).nm,d_,n_,tissue(1).nm, ...
            sim(1),sim(2),sim(3),n_dz,sim(4),sim(5),size(tissue,2));
            for layer = 1:size(tissue,2)
                fprintf(fileID,['%f	%f	%f	%f	%f    	# layer %d\n'], ...
                    tissue(layer).n, ...
                    tissue(layer).mua, ...
                    tissue(layer).mus, ...
                    tissue(layer).g, ...
                    tissue(layer).d, ...
                    layer);
            end           
            fprintf(fileID,['1                        	# n for medium below\n\n']);
    end
end


   


