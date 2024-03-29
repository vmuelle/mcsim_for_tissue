function getmcml(name)
% function getmcml(name)
% 	reads the MCML output file <name> and returns the global values.
% 	If PRINTON ==1 or PLOTON ==1, output is turned on.
% USES
%	makec2f.m which creates a colormap.
% by Steven L. Jacques, Jan. 2008
% Oregon Health & Sciences University, Portland, OR, USA

global A Al Az Azr Fzr Fz Na Nc Nlayers 
global Nr Nz Rd Ra Rr Rra Rsp T Ta Td Tr Tra d dr dz 
global g mua mus n nabove nbelow r z
global rm Azrm Fzrm
global PLOTON PRINTON

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
if PRINTON
    disp(sprintf('dr = %0.4f',dr))
    disp(sprintf('dz = %0.4f',dz))
    disp(sprintf('Nr = %d',Nr))
    disp(sprintf('Nz = %d',Nz))
    disp(sprintf('Na = %d',Na))
    disp(sprintf('Nlayer = %d', Nlayers))
end
line = fgetl(fid);
line = fgetl(fid);
nabove = sscanf(line, '%f');
if PRINTON
    disp(sprintf('nabove = %0.3f', nabove))
end

n = zeros(Nlayers,1); mua = n; mus = n; g = n; d = n;
if PRINTON
    disp(sprintf('\t#:\tn   \tmua  \tmus  \tg    \td'))
end
for i=1:Nlayers
	line = fgetl(fid);
	u = sscanf(line, '%f %f %f %f %f');
	n(i)   = u(1);
	mua(i) = u(2);
	mus(i) = u(3);
	g(i)   = u(4);
	d(i)   = u(5);
	if PRINTON
        disp(sprintf('\t%d:\t%0.2f\t%0.2f\t%0.1f\t%0.3f\t%0.4f', i, n(i), mua(i), mus(i), g(i), d(i) ))
    end
end
line = fgetl(fid);
nbelow = sscanf(line,'%f');
if PRINTON
    disp(sprintf('nbelow = %0.3f', nbelow))
end
line = fgetl(fid);

%%%%
% read RAT ---> Rsp, Rd, A, T
%%%%
line = fgetl(fid); 
if PRINTON
    disp(line)
end
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
if PRINTON
    disp(line)
end
for i=1:Nlayers
	line = fgetl(fid);
	Al(i) = sscanf(line,'%f');
	if PRINTON
        disp(sprintf('\tAl(%d) = %f',i,  Rsp))
    end
end

%%%
% Read A_z --> Az
%%%%
line = fgetl(fid); 
line = fgetl(fid); 
if PRINTON
    disp(line)
end
Az = fscanf(fid, '%f');


%%%
% read Rd_r ---> Rr[Nr]
%%%%
line = fgetl(fid); 
if PRINTON
    disp(line)
end
Rr = fscanf(fid,'%f');

%%%
% read Rd_a ---> Ra[Na]
%%%%
line = fgetl(fid); disp(line)
for i=1:Na
	line = fgetl(fid);
	Ra(i,1) = sscanf(line,'%f');
end

%%%
% read Tt_r ---> Tr[Nr]
%%%%
line = fgetl(fid); 
line = fgetl(fid); 
if PRINTON
    disp(line)
end
Tr = fscanf(fid, '%f');

%%%
% read Tt_a ---> Ta[Na]
%%%%
line = fgetl(fid); 
if PRINTON
    disp(line)
end
for i=1:Na
	line = fgetl(fid);
	Ta(i,1) = sscanf(line,'%f');
end

%%%
% read A_rz ---> Azr[Nz,Nr]
%%%%
line = fgetl(fid); 
line = fgetl(fid); 
if PRINTON
    disp(line)
end
for i=1:5; line = fgetl(fid); end
u     = fscanf(fid,'%f');
Azr   = reshape(u, [Nz Nr]);

%%%
% read Rd_ra ---> Rra[Nr, Na]
%%%%
line = fgetl(fid); 
if PRINTON
    disp(line)
end
for i=1:5; line = fgetl(fid); end
u = fscanf(fid, '%f');
Rra = reshape(u, [Nr Na]);

%%%
% read Tt_ra ---> Tra[Nr, Na]
%%%%
line = fgetl(fid); 
if PRINTON
    disp(line)
end
for i=1:5; line = fgetl(fid); end
u = fscanf(fid, '%f');
Tra = reshape(u, [Nr Na]);

%%%
% done
%%%
fclose(fid);
clear fid


%%%%
% r, z
%%%
z     = ((1:Nz)' - 0.5)*dz;
r     = ((1:Nr)' - 0.5)*dr;

%%%
% Azrm
%%%
Azrm  = zeros(Nz, 2*Nr-1);
Azrm(:,Nr:2*Nr-1) = Azr;
Azrm(:,1:Nr) = Azr(:,Nr:-1:1);
rm = [((-Nr+1:-1) - 0.5)*dr (r'- dr)]' + dr/2;
% rm was checked by plot(diff(rm)) and imagesc(rm,z,log10(Azrm))

% find ma(iz)
j = 1;
D = d(j);
for iz=1:Nz
	if z(iz) < D
		ma(iz) = mua(j);
	else
		j = j+1;
		D = D + d(j);
		ma(iz) = mua(j);
	end
end

Fzrm = Azrm*0; % set size of Fzrm
Fz = zeros(Nz,1);
for iz=1:Nz
	Fzrm(iz,:) = Azrm(iz,:)/ma(iz);
    Fzr(iz,:) = Azr(iz,:)/ma(iz);
	Fz(iz) = Az(iz)/ma(iz);
end

if PLOTON % 0 = plotting turned OFF
%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%
%%%%%%
% Plot Azr, Fzr
%%%
figure(1);clf
figure(2);clf
figure(3);clf

set(figure(1),'position',[10    15   700   800])
set(figure(2),'position',[735   476   600   360])
set(figure(3),'position',[735    15   600   360])

sz = 12;
u = (2:length(rm)-2);
v = (1:Nz-1);

figure(1)
subplot(2,1,1)
imagesc(rm(u),z(v),log10(Azrm(v,u)));
set(gca,'fontsize',sz)
xlabel('r [cm]')
ylabel('z [cm]')
title(sprintf('log_1_0( Azr [J/cm^3] ),       Rd = %0.5f', Rd))
colorbar
colormap(makec2f)
set(colorbar,'fontsize',sz)
axis('equal')

subplot(2,1,2)
imagesc(rm(u),z(v),log10(Fzrm(v,u)));
set(gca,'fontsize',sz)
xlabel('r [cm]')
ylabel('z [cm]')
title('log_1_0( Fzr [J/cm^2] )')
colorbar
colormap(makec2f) % <--------------- need makec2f.m
set(colorbar,'fontsize',sz)
axis('equal')

%%%%
% plot Rr
%%%
figure(2)
u = (1:Nr-1);
semilogy(r(u), Rr(u),'r-')
hold on
plot(r(u), Tr(u),'b-')
plot(-r(u), Rr(u),'r-')
hold on
plot(-r(u), Tr(u),'b-')
set(gca,'fontsize',sz)
xlabel('r [cm]')
ylabel('R [cm^-^2]')
title(name)
legend('Rr','Tr')

%%%%
% plot Az
%%%
figure(3)
u = (1:Nz-1);
subplot(2,1,1)
	plot(z(u), Az(u),'ro-')
	set(gca,'fontsize',sz)
	xlabel('r [cm]')
	ylabel('Az [cm^-^1]')
subplot(2,1,2)
	plot(z(u), Fz(u),'ro-')
	set(gca,'fontsize',sz)
	xlabel('r [cm]')
	ylabel('Fz [-]')

end %%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%

clear D ans a b c2 h h2 i ii j line sz u v w w2 ifile
