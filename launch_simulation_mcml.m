function launch_simulation_mcml()
    clear all;
 
    %get the parameters for tissue with different properties
    tissue = Tissue();
    tissue_list_dn = tissue.makeTissueList('diastolic','normal');
    tissue_list_dc = tissue.makeTissueList('diastolic','compressed');
    tissue_list_sn = tissue.makeTissueList('systolic','normal');
    tissue_list_sc = tissue.makeTissueList('systolic','compressed');
    
    
    %set parameters for simulation 
    photons = 10000000;             %no of photons
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
    launch_simulation_single(filename);
    
    filename = 'data_files/inputs/coefficients_d_c.mci';
    fileID = fopen(filename,'w');
    print_to_file(tissue_list_dc,fileID,'d','c',sim);
    fclose(fileID);
    launch_simulation_single(filename);
    
    filename = 'data_files/inputs/coefficients_s_n.mci';
    fileID = fopen(filename,'w');
    print_to_file(tissue_list_sn,fileID,'s','n',sim);
    fclose(fileID);
    launch_simulation_single(filename);
    
    filename = 'data_files/inputs/coefficients_s_c.mci';
    fileID = fopen(filename,'w');
    print_to_file(tissue_list_sc,fileID,'s','c',sim);
    fclose(fileID);
    launch_simulation_single(filename);
end

function launch_simulation_single(filename)
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
