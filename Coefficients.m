%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculation for the coefficients of absorption and scattering for 
% different skin tissues to input into mcml 
%
% Autor:         Viktor Müller
% Datum:         21.01.22
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Coefficients < handle
   
    properties(Access = private)
        %Given parameters for skin layers normal
        n_n;
        d_d_n;
        c_b_n;
        c_w_n;
        c_f_n;
        v_d_n;
        sc_n;
        p_n;

        %Given parameters for skin layers compressed
        n_c;
        d_d_c;
        c_b_c;
        c_w_c;
        c_f_c;
        v_d_c;
        sc_c;
        p_c;

        %Parameters for development
        num_freq;
        
        %Fixed parameters from literature
        c_w_0;
        spO2;
        svO2;
        p;
       
        ratio_a_nv_n; 
        ratio_a_nv_c;
        absorption_water; 
        absorption_hb;  
        absorption_hbO2; 
        absorption_fat; 
        sc_b; 
        
        %Matrices to store the results and intermediate results
        %used for absorption
        absorption_complete_d_n;         %diastolic normal
        absorption_complete_s_n;         %systolic normal 
        absorption_complete_d_c;         %disatolic compressed skin
        absorption_complete_s_c;         %systolic compressed skin
        absorption_base_n;  
        absorption_base_c; 
        c_b_d_n_corrected; 
        c_b_d_c_corrected; 
        c_b_s_n_corrected; 
        c_b_s_c_corrected; 
        f_a_n;  
        f_v_n; 
        f_a_c;  
        f_v_c;  
        d_s_n; 
        d_s_c;
        f_a_s_n; 
        f_a_s_c; 
        
        %used for scattering
        scattering_n; 
        scattering_c; 
        
        %used for anisotrophy
        anisotrophy;
    end
    properties(Access = public)
        lambda; 
        complete_coefficients_normal_systolic;
        complete_coefficients_compressed_systolic;
        complete_coefficients_normal_diastolic;
        complete_coefficients_compressed_diastolic;
    end

    methods
        function obj = Coefficients(normal_skin,compressed_skin,fixed_parameters,freqency_dependent_parameters)
            %save all parameters into class variables
            obj.n_n = normal_skin(1,:);
            obj.d_d_n = normal_skin(2,:);
            obj.c_b_n = normal_skin(3,:);
            obj.c_w_n = normal_skin(4,:);
            obj.c_f_n = normal_skin(5,:);
            obj.v_d_n = normal_skin(6,:);
            obj.sc_n = normal_skin(7,:);
            obj.ratio_a_nv_n = normal_skin(8,:);
            obj.p_n = normal_skin(9,:);
            obj.n_c = compressed_skin(1,:);
            obj.d_d_c = compressed_skin(2,:);
            obj.c_b_c = compressed_skin(3,:);
            obj.c_w_c = compressed_skin(4,:);
            obj.c_f_c = compressed_skin(5,:);
            obj.v_d_c = compressed_skin(6,:);
            obj.sc_c = compressed_skin(7,:);
            obj.ratio_a_nv_c = compressed_skin(8,:);
            obj.p_c = compressed_skin(9,:);
            obj.c_w_0 = fixed_parameters(1);
            obj.spO2 = fixed_parameters(2);
            obj.svO2 = fixed_parameters(3);
            obj.p = fixed_parameters(4);
            obj.anisotrophy= fixed_parameters(5);
            obj.lambda = freqency_dependent_parameters(1,:);
            obj.absorption_water = freqency_dependent_parameters(2,:);
            obj.absorption_hb = freqency_dependent_parameters(3,:);
            obj.absorption_hbO2 = freqency_dependent_parameters(4,:);
            obj.absorption_fat = freqency_dependent_parameters(5,:);
            obj.sc_b = freqency_dependent_parameters(6,:);
            obj.num_freq = size(obj.lambda,2);

            %initialize all unknown classvariables with zero
            obj.absorption_complete_d_n = zeros(6,obj.num_freq);        %diastolic normal
            obj.absorption_complete_s_n = zeros(6,obj.num_freq);        %systolic normal 
            obj.absorption_complete_d_c = zeros(5,obj.num_freq);        %disatolic compressed skin
            obj.absorption_complete_s_c = zeros(5,obj.num_freq);        %systolic compressed skin
            obj.absorption_base_n = zeros(6,obj.num_freq);
            obj.absorption_base_c = zeros(5,obj.num_freq);
            obj.c_b_d_n_corrected = zeros(6,obj.num_freq);
            obj.c_b_d_c_corrected = zeros(5,obj.num_freq);
            obj.c_b_s_n_corrected = zeros(6,obj.num_freq);
            obj.c_b_s_c_corrected = zeros(5,obj.num_freq);
            obj.f_a_n = zeros (6,obj.num_freq);
            obj.f_v_n = zeros (6,obj.num_freq);
            obj.f_a_c = zeros (6,obj.num_freq);
            obj.f_v_c = zeros (6,obj.num_freq);
            obj.d_s_n = zeros (6,obj.num_freq);
            obj.d_s_n(1,:) = obj.d_d_n(1);
            obj.d_s_c = zeros (5,obj.num_freq);
            obj.d_s_c(1,:) = obj.d_d_c(1);
            obj.f_a_s_n = zeros(6,obj.num_freq);
            obj.f_a_s_c = zeros(5,obj.num_freq);
            obj.complete_coefficients_normal_systolic = zeros(6,obj.num_freq,5);
            obj.complete_coefficients_compressed_systolic = zeros(5,obj.num_freq,5);
            obj.complete_coefficients_normal_diastolic = zeros(6,obj.num_freq,5);
            obj.complete_coefficients_compressed_diastolic = zeros(5,obj.num_freq,5);
        
            %used for scattering
            obj.scattering_n = zeros(6,obj.num_freq);
            obj.scattering_c = zeros(5,obj.num_freq);


            obj.get_absorption_from_extrinction();
            obj.blood_concentration_effective_diastolic();
            obj.blood_concentration_effective_systolic();
            obj.absorption_first_layer();
            obj.absorption_base();
            obj.absorption_without_blood();
            obj.absorption_systolic_state();
            obj.absorption_diastolic_state();
            obj.scattering();
            obj.make_final_coefficient_arrays();

        end
    
        function get_absorption_from_extrinction(self)
        %transform extrinction coefficients into absorption coefficients
            for freq = 1:self.num_freq
                self.absorption_hb(freq) = 2.303 * self.absorption_hb(freq) * 150 / 64500;
                self.absorption_hbO2(freq) = 2.303 * self.absorption_hbO2(freq) * 150 / 64500;
            end
        end
        
       
        function blood_concentration_effective_diastolic(self)
        %determine effective blood concentraition with different parameters for
        %arterial and venous blood in diastolic state
            for freq = 1:self.num_freq
                for layer = 2:6
                    self.f_a_n(layer,freq) = self.ratio_a_nv_n(layer) * self.c_b_n(layer) * self.arterial_venous_c_b(((1-self.spO2)*self.absorption_hb(freq))+(self.spO2*self.absorption_hbO2(freq)),self.v_d_n(layer));
                    self.f_v_n(layer,freq) = (1-self.ratio_a_nv_n(layer)) * self.c_b_n(layer) * self.arterial_venous_c_b(((1-self.svO2)*self.absorption_hb(freq))+(self.svO2*self.absorption_hbO2(freq)),self.v_d_n(layer));
                    self.c_b_d_n_corrected(layer,freq) = self.f_a_n(layer,freq) + self.f_v_n(layer,freq);
                end
                
                for layer = 2:5 
                    self.f_a_c(layer,freq) = self.ratio_a_nv_c(layer) * self.c_b_c(layer) * self.arterial_venous_c_b(((1-self.spO2)*self.absorption_hb(freq))+(self.spO2*self.absorption_hbO2(freq)),self.v_d_c(layer)); 
                    self.f_v_c(layer,freq) = (1-self.ratio_a_nv_c(layer)) * self.c_b_c(layer) * self.arterial_venous_c_b(((1-self.svO2)*self.absorption_hb(freq))+(self.svO2*self.absorption_hbO2(freq)),self.v_d_c(layer)); 
                    self.c_b_d_c_corrected(layer,freq) = self.f_a_c(layer,freq) + self.f_v_c(layer,freq);
                end
            end
        end
        
       
        function blood_concentration_effective_systolic(self)
        %determine effective blood concentraition with different parameters for
        %arterial and venous blood in systolic state
            for freq = 1:self.num_freq
                for layer = 2:6
                    self.f_a_s_n(layer,freq) = self.f_a_n(layer,freq) + self.p_n(layer)*self.arterial_venous_c_b(((1-self.spO2)*self.absorption_hb(freq))+(self.spO2*self.absorption_hbO2(freq)),self.v_d_n(layer));
                    self.c_b_s_n_corrected(layer,freq) = self.c_b_d_n_corrected(layer,freq)+self.p_n(layer)*self.arterial_venous_c_b(((1-self.spO2)*self.absorption_hb(freq))+(self.spO2*self.absorption_hbO2(freq)),self.v_d_n(layer));
                    self.d_s_n(layer,freq) = self.d_d_n(layer)*(1+self.c_b_s_n_corrected(layer,freq)*self.p_n(layer));
                end
                
                for layer = 2:5 
                    self.f_a_s_c(layer,freq) = self.f_a_c(layer,freq) + self.p_c(layer)*self.arterial_venous_c_b(((1-self.spO2)*self.absorption_hb(freq))+(self.spO2*self.absorption_hbO2(freq)),self.v_d_c(layer));
                    self.c_b_s_c_corrected(layer,freq) = self.c_b_d_c_corrected(layer,freq)+self.p_c(layer)*self.arterial_venous_c_b(((1-self.spO2)*self.absorption_hb(freq))+(self.spO2*self.absorption_hbO2(freq)),self.v_d_c(layer));
                    self.d_s_c(layer,freq) = self.d_d_c(layer)*(1+self.c_b_s_c_corrected(layer,freq)*self.p_c(layer));
                end
            end
        end
        
        
       
        function absorption_first_layer(self)
        %calculate absorption coefficients for first layer (epidermis)
            for freq = 1:self.num_freq
                self.absorption_complete_d_n(1,freq) = (self.c_w_n(1)*self.absorption_water(freq)) + (1-self.c_w_n(1)) * 0.5*(0.244+(85.3*exp((-self.lambda(freq)-154)/66.2)));
                self.absorption_complete_s_n(1,freq) = self.absorption_complete_d_n(1,freq);
            
                self.absorption_complete_d_c(1,freq) = (self.c_w_c(1)*self.absorption_water(freq)) + (1-self.c_w_c(1)) * 0.5*(0.244+(85.3*exp((-self.lambda(freq)-154)/66.2)));
                self.absorption_complete_s_c(1,freq) = self.absorption_complete_d_c(1,freq);
            end 
        end
        
       
        function absorption_base(self)
        %calculate base absorption coefficients (normal and compressed)
            for freq = 1:self.num_freq
                for layer = 2:5
                    self.absorption_base_n(layer,freq) = 0.5 * (self.c_w_n(layer)/self.c_w_0)* (0.244 + 16.82* exp((-self.lambda(freq)-400)/80.5));
                end
                self.absorption_base_n(6,freq) = 0.25 * (self.c_w_n(layer)/self.c_w_0)* (0.244 + 16.82* exp((-self.lambda(freq)-400)/80.5));
                
                for layer = 2:4
                    self.absorption_base_c(layer,freq) = 0.5 * (self.c_w_n(layer)/self.c_w_0)* (0.244 + 16.82* exp((-self.lambda(freq)-400)/80.5));
                end
                self.absorption_base_c(5,freq) = 0.25 * (self.c_w_n(layer)/self.c_w_0)* (0.244 + 16.82* exp((-self.lambda(freq)-400)/80.5));
            end
        end
        
       
        function absorption_without_blood(self)
        %calculate absorption coefficients for all other layers WITHOUT blood (normal and
        %compressed)
            for freq = 1:self.num_freq
                for layer = 2:6
                    self.absorption_complete_d_n(layer,freq) = self.c_f_n(layer)*self.absorption_fat(freq) + (1-self.c_f_n(layer))* self.c_w_n(layer) * self.absorption_water(freq) + (1-self.c_f_n(layer))*(1-self.c_w_n(layer))*self.absorption_base_n(layer,freq);
                    self.absorption_complete_s_n(layer,freq) = self.c_f_n(layer)*self.absorption_fat(freq) + (1-self.c_f_n(layer))* self.c_w_n(layer) * self.absorption_water(freq) + (1-self.c_f_n(layer))*(1-self.c_w_n(layer))*self.absorption_base_n(layer,freq);
                end
                
                for layer = 2:5
                    self.absorption_complete_d_c(layer,freq) = self.c_f_c(layer)*self.absorption_fat(freq) + (1-self.c_f_c(layer))* self.c_w_c(layer) * self.absorption_water(freq) + (1-self.c_f_c(layer))*(1-self.c_w_c(layer))*self.absorption_base_c(layer,freq);
                    self.absorption_complete_s_c(layer,freq) = self.c_f_c(layer)*self.absorption_fat(freq) + (1-self.c_f_c(layer))* self.c_w_c(layer) * self.absorption_water(freq) + (1-self.c_f_c(layer))*(1-self.c_w_c(layer))*self.absorption_base_c(layer,freq);
                end
            end
        end
        
        %calculate absorption coefficients for all other layers WITH blood (normal and
        %compressed) in diastolic state
        function absorption_diastolic_state(self)
            for freq = 1:self.num_freq
                for layer = 2:6
                    self.absorption_complete_d_n(layer,freq) = self.f_a_n(layer,freq)*((1-self.spO2)*self.absorption_hb(freq)+self.spO2*self.absorption_hbO2(freq)) ...
                        + self.f_v_n(layer,freq)*((1-self.svO2)*self.absorption_hb(freq)+self.svO2*self.absorption_hbO2(freq)) ...
                        + self.c_b_d_n_corrected(layer,freq)*self.absorption_water(freq)+(1-self.c_b_n(layer))*self.absorption_complete_d_n(layer,freq);
                end
                
                for layer = 2:5
                    self.absorption_complete_d_c(layer,freq) = self.f_a_c(layer,freq)*((1-self.spO2)*self.absorption_hb(freq)+self.spO2*self.absorption_hbO2(freq))... 
                        + self.f_v_c(layer,freq)*((1-self.svO2)*self.absorption_hb(freq)+self.svO2*self.absorption_hbO2(freq))...
                        + self.c_b_d_c_corrected(layer,freq)*self.absorption_water(freq)+(1-self.c_b_c(layer))*self.absorption_complete_d_c(layer,freq);
                end
            end
        end
        
       
        function absorption_systolic_state(self)
        %calculate absorption coefficients for all other layers WITH blood (normal and
        %compressed) in systolic state
            for freq = 1:self.num_freq
                for layer = 2:6
                    self.absorption_complete_s_n(layer,freq) = self.d_d_n(layer)/self.d_s_n(layer,freq) ...
                        * (... 
                        self.f_a_s_n(layer,freq)*((1-self.spO2)*self.absorption_hb(freq)+self.spO2*self.absorption_hbO2(freq)) ...
                        + self.f_v_n(layer,freq)*((1-self.svO2)*self.absorption_hb(freq)+self.svO2*self.absorption_hbO2(freq))...
                        + self.c_b_s_n_corrected(layer,freq)*self.absorption_water(freq)+(1-self.c_b_s_n_corrected(layer,freq))*self.absorption_complete_s_n(layer,freq)...
                        );
                end
                
                for layer = 2:5
                    self.absorption_complete_s_c(layer,freq) = self.d_d_c(layer)/self.d_s_c(layer,freq) ...
                        * (... 
                        self.f_a_s_c(layer,freq)*((1-self.spO2)*self.absorption_hb(freq)+self.spO2*self.absorption_hbO2(freq)) ...
                        + self.f_v_c(layer,freq)*((1-self.svO2)*self.absorption_hb(freq)+self.svO2*self.absorption_hbO2(freq))...
                        + self.c_b_s_c_corrected(layer,freq)*self.absorption_water(freq)+(1-self.c_b_s_c_corrected(layer,freq))*self.absorption_complete_s_c(layer,freq)...
                        );
                end
            end
        end

        
        function scattering(self)
        %Get scattering parameters
            for freq = 1:self.num_freq
                for layer = 1:6
                    self.scattering_n(layer,freq) = self.sc_n(layer)* self.lambda(freq)^(-self.sc_b(freq)); 
                end
                for layer = 1:5
                    self.scattering_c(layer,freq) = self.sc_c(layer)* self.lambda(freq)^(-self.sc_b(freq));
                end
            end
        end
    
        function c_b_x = arterial_venous_c_b(self,absorption,v_d)
        %calculate arterial and venous c
            c_b_x = 1/(1+1.007*(absorption*v_d/2)^1.228);
        end

        function make_final_coefficient_arrays(self)
        %write coefficients into a sturcure, that has every parameter
        %needed for the MCML input in a 2D array and for ever frequency in
        %the third dimension
            for freq = 1:self.num_freq
                for layer = 1:6
                    self.complete_coefficients_normal_diastolic(layer,freq,:) = ...
                        [self.n_n(layer) self.absorption_complete_d_n(layer,freq) self.scattering_n(layer,freq) self.anisotrophy self.d_d_n(layer)];
                    self.complete_coefficients_normal_systolic(layer,freq,:) = ...
                        [self.n_n(layer) self.absorption_complete_s_n(layer,freq) self.scattering_n(layer,freq) self.anisotrophy self.d_s_n(layer,freq)];
                end
                for layer = 1:5
                    self.complete_coefficients_compressed_diastolic(layer,freq,:) = ...
                        [self.n_c(layer) self.absorption_complete_d_c(layer,freq) self.scattering_c(layer,freq) self.anisotrophy self.d_d_c(layer)];
                    self.complete_coefficients_compressed_systolic(layer,freq,:) = ...
                        [self.n_c(layer) self.absorption_complete_s_c(layer,freq) self.scattering_c(layer,freq) self.anisotrophy self.d_s_c(layer,freq)];
                end
            end
        end
    end
end