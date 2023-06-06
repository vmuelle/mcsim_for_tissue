    
function vismcml1photon(n)
    tissue = Tissue();
    p = 0.3;
    tissue_list_dn = tissue.makeTissueList('diastolic','normal',p);
    photons = 1;             %no of photons
    dz = 0.001;                   %separation between grid lines (cm) in z direction                                                                                                                                                                 
    dr = 0.001;                   %separation between grid lines (cm) in r direction
    n_dr = 1000;                  %no of grid elements in r direction        
    n_da = 90;                    %no of parts of an 90 degree angle, that is spanned between the photon exiting direction and the surface normal
    sim = [photons,dz,dr,n_dr,n_da];
    
    figure
    %print mcml input files and launch simulation
    filename = 'data_files/inputs/coefficients_d_n_single.mci';
    fileID = fopen(filename,'w');
    print_to_file(tissue_list_dn,fileID,'d','n',sim);
    fclose(fileID);
    for i = 1:n
        launch_simulation_single(filename);
        getmcml('data_files/outputs/single_photon_d_n_470.mco');
    end
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
            'data_files/outputs/single_photon_%c_%c_%d.mco 	A            	# output file name, ASCII.\n'...
            '%d	                  	    # No. of photons\n'...
            '%f	%f               	# dz, dr [cm]\n'...
            '%d	%d	%d                	# No. of dz, dr, da.\n' ...
            '\n'...
            '%d                        	# Number of layers\n'...
            '#n	mua	mus	g	d         	# One line for each layer\n'...
            '1.0                         	# n for medium above\n'], ...
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
            fprintf(fileID,['1.0                        	# n for medium below\n\n']);
    end
end






function getmcml(name)
% function getmcml(name)
% 	reads the MCML output file <name> and returns the global values.
% 	If PRINTON ==1 or PLOTON ==1, output is turned on.
% USES
%	makec2f.m which creates a colormap.
% by Steven L. Jacques, Jan. 2008
% Oregon Health & Sciences University, Portland, OR, USA

PLOTON = 1;
PRINTON = 0;

fid = fopen(name,'r');

% PLOTON = 0; % 0 = off, 1 = on.  Controls plotting option.
% PRINTON = 0;

%%%
% read dr, dz, Nz, Nr, Na, Nlayer
%%%%
for i=1:15; line = fgetl(fid); disp(line); end
u = sscanf(line, '%f %f');
dz = u(1);
dr = u(2);
line = fgetl(fid);
u = sscanf(line, '%d %d %d');
Nz = u(1);
Nr = u(2);
Na = u(3);
line = fgetl(fid);
line = fgetl(fid);
Nlayers = sscanf(line, '%d');

line = fgetl(fid);
line = fgetl(fid);
nabove = sscanf(line, '%f');

n = zeros(Nlayers,1); mua = n; mus = n; g = n; d = n;

for i=1:Nlayers
	line = fgetl(fid);
	u = sscanf(line, '%f %f %f %f %f');
	n(i)   = u(1);
	mua(i) = u(2);
	mus(i) = u(3);
	g(i)   = u(4);
	d(i)   = u(5);
end
line = fgetl(fid);
nbelow = sscanf(line,'%f');
line = fgetl(fid);

%%%%
% read RAT ---> Rsp, Rd, A, T
%%%%
line = fgetl(fid); 
line = fgetl(fid);
Rsp  = sscanf(line,'%f'); disp(sprintf('\tRsp = %f', Rsp))
line = fgetl(fid);
Rd    = sscanf(line,'%f'); disp(sprintf('\tRd = %f', Rd))
line = fgetl(fid);
A    = sscanf(line,'%f'); disp(sprintf('\tA = %f', A))
line = fgetl(fid);
Td   = sscanf(line,'%f'); disp(sprintf('\tT = %f', Td))

%%%
% Read A_l --> Al
%%%%
line = fgetl(fid);
line = fgetl(fid); 
for i=1:Nlayers
	line = fgetl(fid);
	Al(i) = sscanf(line,'%f');
end

% Read A_z --> Az
line = fgetl(fid); 
line = fgetl(fid); 
Az = fscanf(fid, '%f');


% read Rd_r ---> Rr[Nr]
line = fgetl(fid); 
Rr = fscanf(fid,'%f');

% read Rd_a ---> Ra[Na]
line = fgetl(fid); disp(line)
for i=1:Na
	line = fgetl(fid);
	Ra(i,1) = sscanf(line,'%f');
end

% read Tt_r ---> Tr[Nr]
line = fgetl(fid); 
line = fgetl(fid); 
Tr = fscanf(fid, '%f');

% read Tt_a ---> Ta[Na]
line = fgetl(fid); 
for i=1:Na
	line = fgetl(fid);
	Ta(i,1) = sscanf(line,'%f');
end

% read A_rz ---> Azr[Nz,Nr]
line = fgetl(fid); 
line = fgetl(fid); 
for i=1:5; line = fgetl(fid); end
u     = fscanf(fid,'%f');
Azr   = reshape(u, [Nz Nr]);

% done
fclose(fid);
clear fid

% r, z
z     = ((1:Nz)' - 0.5)*dz;
r     = ((1:Nr)' - 0.5)*dr;


if PLOTON % 0 = plotting turned OFF
%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%
%%%%%%
% Plot Azr
%%%

u = (2:length(r)-2);
v = (1:Nz-1);

%toplot = Azrm(v,u) > 0;
sort_azr = sort(reshape(Azr,1,[]),'descend');
hold on
j = 1;
for i = 1:length(sort_azr)
    if sort_azr(i) ~= 0
        [photon_index_row(j),photon_index_column(j)] = find(Azr==sort_azr(i));
        j = j+1;
    end
end
plot(photon_index_column.*dr,photon_index_row.*dz,'-o');

set(gca, 'YDir','reverse')
xlabel('r [cm]')
xlim([0 1])
ylabel('z [cm]')
ylim([0 0.8])
title('Tracing some Photons')

end %%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%
end

