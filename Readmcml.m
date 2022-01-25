classdef Readmcml < handle
    properties(Access = private)
        A; Al; Azr; Fzr; Fz; Na; Nc; Nlayers; 
        Nr; Nz; Rd; Ra; Rr; Rra; Rsp; T; Ta; Td; Tr; Tra; dr; 
        g; mua; mus; n; nabove; nbelow; r; z;
        rm; Azrm; Fzrm;
    end
    properties(Access = public)
        Az;d;dz;
    end
    methods
        function obj = Readmcml(name)
            obj.readmcml(name);
        end
        function readmcml(self,name)
        % function getmcml(name)
        % 	reads the MCML output file <name> and stores them into class
        % 	properties 
        
        fid = fopen(name,'r');
        
        
        %%%
        % read dr, dz, Nz, Nr, Na, Nlayer
        %%%%
        for i=1:15; 
            line = fgetl(fid);    
        end
        u = sscanf(line, '%f %f');
        self.dz = u(1);
        self.dr = u(2);
        line = fgetl(fid);
        u = sscanf(line, '%d %d %d');
        self.Nz = u(1);
        self.Nr = u(2);
        self.Na = u(3);
        line = fgetl(fid);
        line = fgetl(fid);
        self.Nlayers = sscanf(line, '%d');
        
        line = fgetl(fid);
        line = fgetl(fid);
        self.nabove = sscanf(line, '%f');
        
        
        self.n = zeros(self.Nlayers,1); 
        self.mua = self.n;
        self.mus = self.n; 
        self.g = self.n; 
        self.d = self.n;
        
        for i=1:self.Nlayers
	        line = fgetl(fid);
	        u = sscanf(line, '%f %f %f %f %f');
	        self.n(i)   = u(1);
	        self.mua(i) = u(2);
	        self.mus(i) = u(3);
	        self.g(i)   = u(4);
	        self.d(i)   = u(5);
	        
        end
        line = fgetl(fid);
        self.nbelow = sscanf(line,'%f');
        
        line = fgetl(fid);
        
        %%%%
        % read RAT ---> Rsp, Rd, A, T
        %%%%
        line = fgetl(fid); 
        line = fgetl(fid);
        self.Rsp  = sscanf(line,'%f'); 
        line = fgetl(fid);
        self.Rd    = sscanf(line,'%f'); 
        line = fgetl(fid);
        self.A    = sscanf(line,'%f'); 
        line = fgetl(fid);
        self.Td   = sscanf(line,'%f'); 
        
        %%%
        % Read A_l --> Al
        %%%%
        line = fgetl(fid);
        line = fgetl(fid); 
        for i=1:self.Nlayers
	        line = fgetl(fid);
	        self.Al(i) = sscanf(line,'%f');
        end
        
        %%%
        % Read A_z --> Az
        %%%%
        line = fgetl(fid); 
        line = fgetl(fid); 
        self.Az = fscanf(fid, '%f');
        
        
        %%%
        % read Rd_r ---> Rr[Nr]
        %%%%
        line = fgetl(fid); 
        self.Rr = fscanf(fid,'%f');
        
        %%%
        % read Rd_a ---> Ra[Na]
        %%%%
        line = fgetl(fid);
        for i=1:self.Na
	        line = fgetl(fid);
	        self.Ra(i,1) = sscanf(line,'%f');
        end
        
        %%%
        % read Tt_r ---> Tr[Nr]
        %%%%
        line = fgetl(fid); 
        line = fgetl(fid); 
        self.Tr = fscanf(fid, '%f');
        
        %%%
        % read Tt_a ---> Ta[Na]
        %%%%
        line = fgetl(fid); 
        for i=1:self.Na
	        line = fgetl(fid);
	        self.Ta(i,1) = sscanf(line,'%f');
        end
        
        %%%
        % read A_rz ---> Azr[Nz,Nr]
        %%%%
        line = fgetl(fid); 
        line = fgetl(fid); 
        for i=1:5; line = fgetl(fid); end
        u     = fscanf(fid,'%f');
        self.Azr   = reshape(u, [self.Nz self.Nr]);
        
        %%%
        % read Rd_ra ---> Rra[Nr, Na]
        %%%%
        line = fgetl(fid); 
        for i=1:5; line = fgetl(fid); end
        u = fscanf(fid, '%f');
        self.Rra = reshape(u, [self.Nr self.Na]);
        
        %%%
        % read Tt_ra ---> Tra[Nr, Na]
        %%%%
        line = fgetl(fid); 
        for i=1:5; line = fgetl(fid); end
        u = fscanf(fid, '%f');
        self.Tra = reshape(u, [self.Nr self.Na]);
        
        %%%
        % done
        %%%
        fclose(fid);
        clear fid
        
        
        %%%%
        % r, z
        %%%
        self.z     = ((1:self.Nz)' - 0.5)*self.dz;
        self.r     = ((1:self.Nr)' - 0.5)*self.dr;
        
        %%%
        % Azrm
        %%%
        self.Azrm  = zeros(self.Nz, 2*self.Nr-1);
        self.Azrm(:,self.Nr:2*self.Nr-1) = self.Azr;
        self.Azrm(:,1:self.Nr) = self.Azr(:,self.Nr:-1:1);
        self.rm = [((-self.Nr+1:-1) - 0.5)*self.dr (self.r'- self.dr)]' + self.dr/2;
        % rm was checked by plot(diff(rm)) and imagesc(rm,z,log10(Azrm))
        
        % find ma(iz)
        j = 1;
        D = self.d(j);
        for iz=1:self.Nz
	        if self.z(iz) < D
		        ma(iz) = self.mua(j);
	        else
		        j = j+1;
		        D = D + self.d(j);
		        ma(iz) = self.mua(j);
	        end
        end
        
        self.Fzrm = self.Azrm*0; % set size of Fzrm
        self.Fz = zeros(self.Nz,1);
        for iz=1:self.Nz
	        self.Fzrm(iz,:) = self.Azrm(iz,:)/ma(iz);
            self.Fzr(iz,:) = self.Azr(iz,:)/ma(iz);
	        self.Fz(iz) = self.Az(iz)/ma(iz);
        end
    
    
        clear D ans a b c2 h h2 i ii j line sz u v w w2 ifile
        end
    end
end
