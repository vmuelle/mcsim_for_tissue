%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculation for the coefficients of absorption and scattering for 
% different skin tissues to input into mcml and execution by mcml
%
% Autor:         Viktor Müller
% Datum:         09.04.22
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Tissue

    properties(Access = private)
        c_w_0; spO2; svO2; p; anisotrophy; sc_b;
        absorption_hbO2; absorption_hb; absorption_water; absorption_fat; c_b; c_w; c_f; v_d; sc; ratio_a_nv; nm; d;
    end
    methods(Access = public)

        function tissueList = makeTissueList(self,d_,n_)
        %function tissue = makeTissueList(d_,n_)
        %  Returns a list of tissue structures with wavelengths from 400 to 1000 nm
        %  with 50 nm steps (so it's 13 wavelenths) with properties d_: 'diastolic' or 'systolic and n_: 'normal' or
        %  'compressed'
            j = 1;
            %nmLIB = [400 450 500 550 600 650 700 725];
            nmLIB = [470 530 660 770 810 940 1020 1050];
            %for i = 400:50:1000  %for i = 400:50:1000 
            for i = nmLIB
                tissue = self.makeTissue(i,d_,n_);
                tissue(1).nm = i;
                tissueList(j,:) = tissue;
                j=j+1;
            end
        end
    end
    methods(Access = private)

        function tissue = makeTissue(obj,nm_,d_,n_)
        %function tissue = makeTissue(nm_,d_,n_)
        %   Returns the tissue optical properties at the wavelength nm with
        %   properties d_: 'diastolic' or 'systolic and n_: 'normal' or
        %   'compressed'
        %       tissue = [name; n; mua; mus; g; d]';
        %       uses SpectralLIB.mat
              
            obj.nm = nm_;           %wavelenth
            obj.c_w_0 = 0.65;       %coefficient that accounts for background measurement
            obj.spO2 = 0.97;        %aterial oxygen saturation
            obj.svO2 = 0.67;        %venious oxgen saturation
            obj.p = 0.3;            %pulsation increase
            obj.anisotrophy = 0.9;  %anisotrophy
            obj.sc_b = 0.1;         %decaying factor for scattering
            if(obj.nm>580)                                                                          
                obj.sc_b = 0.05;
            end
        
            % Load spectral library
            load spectralLIBv3.mat %load spectralLIB.mat
            %   muadeoxy      701x1              5608  double              
            %   muamel        701x1              5608  double              
            %   muaoxy        701x1              5608  double              
            %   muawater      701x1              5608  double              
            %   musp          701x1              5608  double              
            %   nmLIB         701x1              5608  double              
            obj.absorption_hbO2 = interp1(nmLIB,muaoxy,obj.nm);     %absorption coefficient hb02
            obj.absorption_hb = interp1(nmLIB,muadeoxy,obj.nm);     %absorption coefficient hb
            obj.absorption_water = interp1(nmLIB,muawater,obj.nm);  %absorption coefficient water
            obj.absorption_fat = interp1(nmLIB,muafat,obj.nm);      %absorption coefficient fat
            
            
            
        
            if(strcmp(n_,'normal'))
                %fixed parameters given by moco per layer
                n = [1.33 1.37 1.40 1.40 1.40 1.44];    %reflective index
                obj.d = [0.08 0.015 0.008 0.12 0.05 0.5];   %thickness  
                obj.c_b = [0 0.004 0.02 0.004 0.04 0.03];   %concentration of blood
                obj.c_w = [0.2 0.65 0.65 0.65 0.65 0.05];   %concentration of water
                obj.c_f = [0 0 0 0 0 0.4];                  %concentration of fat
                obj.v_d = [0 10 20 20 40 50];               %vessel diameters
                obj.sc = [15 20 20 20 20 10];               %scattering calibration constant
                obj.ratio_a_nv = [0 0.5 0.5 0.5 0.5 0.5];   %ratio aterial blood (not venious)
                
                tissue(1).name = 'EPIDERMIS';
                tissue(2).name = 'CAPILLARY LOOPS';
                tissue(3).name = 'UPPER PLEXUS';
                tissue(4).name = 'RETICULAR PLEXUS';
                tissue(5).name = 'DEEP PLEXUS';
                tissue(6).name = 'HYPODERMIS';
                
                for i = 1:6
                    tissue(i).n = n(i);
                    tissue(i).d = obj.d(i);
                    tissue(i).g = obj.anisotrophy;
                    %tissue(i).mus = obj.sc(i);
                    tissue(i).mus = obj.scattering(i);
                end
                
                 tissue(1).mua = obj.absorption_first_layer();
        
                for i = 2:6
                    if(strcmp(d_,'diastolic'))
                        tissue(i).mua = obj.absorption_diastolic_state(i);
                    elseif(strcmp(d_,'systolic'))
                        tissue(i).mua = obj.absorption_systolic_state(i);
                    else 
                        error('wrong input paramter for d_')
                    end
                end
        
            elseif(strcmp(n_,'compressed'))
                %fixed paraters, given by moco for compressed skin
                n = [1.33 1.37 1.40 1.40 1.44];         %reflective index
                obj.d = [0.08 0.008 0.004 0.1 0.2];         %thickness
                obj.c_b = [0 0.0012 0.0024 0.024 0.036];    %concentration of blood
                obj.c_w = [0.05 0.15 0.15 0.15 0.35];       %concentration of water
                obj.c_f = [0 0 0 0 0.4];                    %concentration of fat
                obj.v_d = [0 10 20 20 40];                  %vessel diameters
                obj.sc = [15 20 20 20 10];                  %scattering calibration constant
                obj.ratio_a_nv = [0 1 1 1 0.75];            %ratio aterial blood (not venious)
        
                tissue(1).name = 'EPIDERMIS';
                tissue(2).name = 'CAPILLARY LOOPS';
                tissue(3).name = 'UPPER PLEXUS';
                tissue(4).name = 'RETICULAR PLEXUS & DEEP PLEXUS';
                tissue(5).name = 'HYPODERMIS';
        
                for i = 1:5
                    tissue(i).n = n(i);
                    tissue(i).d = obj.d(i);
                    tissue(i).g = obj.anisotrophy;
                    %tissue(i).mus = obj.sc(i);
                    tissue(i).mus = obj.scattering(i);
                end
                
                tissue(1).mua = obj.absorption_first_layer();
        
                for i = 2:5
                    if(strcmp(d_,'diastolic'))
                        tissue(i).mua = obj.absorption_diastolic_state(i);
                    elseif(strcmp(d_,'systolic'))
                        tissue(i).mua = obj.absorption_systolic_state(i);
                    else
                        error('wrong input parameter for d_')
                    end
                end
            else
                error('wrong input parameter for n_')
            end
        end
        
        function [f_a, f_v, c_b_corrected] = blood_concentration_effective_diastolic(obj,layer)
        %function [f_a, f_v, c_b_corrected] = blood_concentration_effective_diastolic(layer)
        % Returns the fraction of aterial/venious blood 'f_a'/'f_v' and the effective
        % blood concentration 'c_b_corrected' for a given layer in diastolc state
        
            f_a = obj.ratio_a_nv(layer) * obj.c_b(layer) * obj.arterial_venous_c_b(((1-obj.spO2)*obj.absorption_hb)+(obj.spO2*obj.absorption_hbO2),obj.v_d(layer));
            f_v = (1-obj.ratio_a_nv(layer)) * obj.c_b(layer) * obj.arterial_venous_c_b(((1-obj.svO2)*obj.absorption_hb)+(obj.svO2*obj.absorption_hbO2),obj.v_d(layer));
            c_b_corrected = f_a + f_v;
        end
        
        
        function [d_s, f_a_s, f_v, c_b_corrected] = blood_concentration_effective_systolic(obj,layer)
        %function [d_s, f_a_s, f_v, c_b_corrected] = blood_concentration_effective_systolic(layer)
        % Returns the corrected thickness 'd_s', the fraction of aterial/venious blood 'f_a_s'/'f_v' and the effective
        % blood concentration 'c_b_corrected' for a given layer in systolc state
        
            [f_a, f_v, c_b_corrected] = obj.blood_concentration_effective_diastolic(layer);
        
            f_a_s = f_a + obj.p*obj.arterial_venous_c_b(((1-obj.spO2)*obj.absorption_hb)+(obj.spO2*obj.absorption_hbO2),obj.v_d(layer));
            c_b_corrected = c_b_corrected+obj.p*obj.arterial_venous_c_b(((1-obj.spO2)*obj.absorption_hb)+(obj.spO2*obj.absorption_hbO2),obj.v_d(layer));
            d_s = obj.d(layer)*(1+c_b_corrected*obj.p);
        end
        
        
        
        function afl =  absorption_first_layer(obj)
        %function afl =  absorption_first_layer()
        % Returns the absorption coefficient 'afl' for first layer (epidermis)
   
            mua_base_epi = 0.5*(0.244+85.3*exp(-(obj.nm-154)/66.2));
            afl = (obj.c_w(1)*obj.absorption_water) + (1-obj.c_w(1)) * mua_base_epi; 
        end
                
               
        function a_base = absorption_base(obj,layer)
        %function a_base = absorption_base(layer)
        % Returns the base absorption 'a_base' for a given 'layer'
              
            factor = 0.5;
            if(layer == length(obj.c_w))
                factor = 0.25;
            end
            a_base = factor * (obj.c_w(layer)/obj.c_w_0)* (0.244 + 16.82* exp(-(obj.nm-400)/80.5));
        end
        
        
        function a_c = absorption_without_blood(obj,layer)
        %function a_c = absorption_without_blood(layer)
        % Returns the absorption coefficient 'a_c' for a given 'layer' WITHOUT blood
        
            a_c = obj.c_f(layer)*obj.absorption_fat + (1-obj.c_f(layer))* obj.c_w(layer) * obj.absorption_water + (1-obj.c_f(layer))*(1-obj.c_w(layer))*obj.absorption_base(layer);
        end
        
        
        function mua = absorption_diastolic_state(obj,layer)
        %function mua = absorption_diastolic_state(layer)
        % Returns the absorption coefficient 'mua' for a given 'layer' WITH
        % blood in diastolic state
        
            [f_a, f_v, c_b_corrected] = obj.blood_concentration_effective_diastolic(layer);
            mua = f_a*((1-obj.spO2)*obj.absorption_hb+obj.spO2*obj.absorption_hbO2) ...
              + f_v*((1-obj.svO2)*obj.absorption_hb+obj.svO2*obj.absorption_hbO2) ...
              + c_b_corrected*obj.absorption_water+(1-obj.c_b(layer))*obj.absorption_without_blood(layer);
        end
        
        
        function mua = absorption_systolic_state(obj,layer)
        %function mua = absorption_systolic_state(layer)
        % Returns the absorption coefficient 'mua' for a given 'layer' WITH blood in systolic state
        
            [d_s, f_a_s, f_v, c_b_corrected] = obj.blood_concentration_effective_systolic(layer);
            mua = obj.d(layer)/d_s ...
                * (... 
                f_a_s*((1-obj.spO2)*obj.absorption_hb+obj.spO2*obj.absorption_hbO2) ...
                + f_v*((1-obj.svO2)*obj.absorption_hb+obj.svO2*obj.absorption_hbO2)...
                + c_b_corrected*obj.absorption_water+(1-c_b_corrected)*obj.absorption_without_blood(layer)...
                );
        end
        
        
        function mus = scattering(obj,layer)
        %function mus = scattering(layer)
        % Returns the scattering coefficient 'mus' for a given 'layer'
        
            mus = obj.sc(layer)* obj.nm^(-obj.sc_b);
        end
        
        function avcb = arterial_venous_c_b(obj,absorption,v_d)
        %function avcb = arterial_venous_c_b(absorption,v_d)
        % Returns a correction factor 'avcb' that compensates the vessel walls to
        % get an apperent blood volume.
        % property 'absorption' hast to be a blood absorption parameter. 
        % property 'v_d' is a vessel diameter
            avcb = 1/(1+1.007*(absorption*v_d/2)^1.228);
        end
    end



end

   


