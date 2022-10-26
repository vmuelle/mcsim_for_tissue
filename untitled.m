%%% Setting up simulation volume %%% 
T = double(ones(10,10,10)); 

T_t = make_tissue(cfg,T,coeff,dz);

radius = 50/2;
center = [200,200,800]/2;
T_tc = make_circle(cfg,T_t,radius,center,3);


function T= make_tissue(T_)
    T = T_;
    lay = 1;
    for iz=1:10 % for every depth z(iz)
        for i = 1:10
            for j = 1:10
                if(lay < cfg.Nt)
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
 
