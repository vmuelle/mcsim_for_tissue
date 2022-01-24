
close all;
clear all;

%Given parameters for skin layers normal
n_n = [1.33 1.37 1.40 1.40 1.40 1.44];
d_d_n = [0.08 0.015 0.008 0.12 0.05 0.5];
c_b_n = [0 0.004 0.02 0.004 0.04 0.03];
c_w_n = [0.2 0.65 0.65 0.65 0.65 0.05];
c_f_n = [0 0 0 0 0 0.4];
v_d_n = [0 10 20 20 40 50];
sc_n = [15 20 20 20 20 10];

ratio_a_nv_n = [0 0.5 0.5 0.5 0.5 0.5];

normal_skin =[n_n; d_d_n; c_b_n; c_w_n; c_f_n; v_d_n; sc_n; ratio_a_nv_n];

%Given parameters for skin layers compressed
n_c = [1.33 1.37 1.40 1.40 1.44];
d_d_c = [0.08 0.008 0.004 0.1 0.2];
c_b_c = [0 0.0012 0.0024 0.024 0.036];
c_w_c = [0.05 0.15 0.15 0.15 0.35];
c_f_c = [0 0 0 0 0.4];
v_d_c = [0 10 20 20 40];
sc_c = [15 20 20 20 10];

ratio_a_nv_c = [0 1 1 1 0.75];

compressed_skin =[n_c; d_d_c; c_b_c; c_w_c; c_f_c; v_d_c; sc_c; ratio_a_nv_c];

%Parameters for development
num_freq= 8;

%Fixed parameters from literature
c_w_0 = 0.65;
spO2 = 0.97;
svO2 = 0.67;
p = 0.1;
anisotrophy = 0.9;
fixed_parameters = [c_w_0; spO2; svO2; p; anisotrophy];

lambda = [400 450 500 550 600 650 700 725];%All following absorption coefficients correspond to these wavelength [nm]
absorption_water = [0.0000663,0.0000922,0.000204,0.000565,0.002224,0.0034,0.00624,0.01489];
absorption_hb = [223296 103292 20862 53412 14677.2 3750.12 1794.28 1244.44]; 
absorption_hbO2 = [266232 62816	20932.8	43016 3200	368 290	364];
absorption_fat = [15 6.375 1.906 0.773 0.466 0.471 0.323 0.415];
sc_b = [0.1 0.1 0.1 0.1 0.05 0.05 0.05 0.05];
freqency_dependent_parameters = [lambda; absorption_water; absorption_hb; absorption_hbO2; absorption_fat; sc_b];


%Parameters for the simulation
photons = 1000;
dz = 0.001;
dr = 0.001;
n_dz = 773;
n_dr = 1000;
n_da = 90;
sim = [photons,dz,dr,n_dz,n_dr,n_da];
%Transform parameters into sufficient values
coefficients = Coefficients(normal_skin,compressed_skin,fixed_parameters,freqency_dependent_parameters);







fileID = fopen('data_files/inputs/coefficients_s_n.mci','w');
print_to_file(coefficients.complete_coefficients_normal_systolic,fileID,coefficients.lambda,'s','n',sim);
fclose(fileID);
fileID = fopen('data_files/inputs/coefficients_d_n.mci','w');
print_to_file(coefficients.complete_coefficients_normal_diastolic,fileID,coefficients.lambda,'d','n',sim);
fclose(fileID);
fileID = fopen('data_files/inputs/coefficients_s_c.mci','w');
print_to_file(coefficients.complete_coefficients_compressed_systolic,fileID,coefficients.lambda,'s','c',sim);
fclose(fileID);
fileID = fopen('data_files/inputs/coefficients_d_c.mci','w');
print_to_file(coefficients.complete_coefficients_compressed_diastolic,fileID,coefficients.lambda,'d','c',sim);
fclose(fileID);

function print_to_file(coeff,fileID,lambda,type_sd,type_cn,sim)
    number_of_runs = size(lambda,2);
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
        '%d                        	# number of runs\n'],type_sd,type_cn,number_of_runs);
    for freq = 1:number_of_runs
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
            freq,lambda(freq),type_sd,type_cn,lambda(freq), ...
            sim(1),sim(2),sim(3),sim(4),sim(5),sim(6),size(coeff,1));
            for layer = 1:size(coeff,1)
                fprintf(fileID,['%f	%f	%f	%f	%f    	# layer %d\n'], ...
                    coeff(layer,freq,1), ...
                    coeff(layer,freq,2), ...
                    coeff(layer,freq,3), ...
                    coeff(layer,freq,4), ...
                    coeff(layer,freq,5), ...
                    layer);
            end           
            fprintf(fileID,[''...
            '1                        	# n for medium below\n\n']);
    end
end


