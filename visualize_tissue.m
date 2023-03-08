function visualize_tissue()
   tissue = Tissue();
   tl = tissue.makeTissueList('diastolic','normal');

   figure
   hold on

   maxmua = 0;
   for i = 1:6
       for j = 1:8
            t = tl(j,i);
            if t.mua > maxmua
                maxmua = t.mua;
            end
            mua(j) = t.mua;  
            lambda(j) = tl(j,1).nm;
       end
       plot(lambda, mua, 'DisplayName', t.name);
   end
   hold off

   xlabel('wavelength [nm]')
   ylabel('mua')
   axis([lambda(1) lambda(end) 0 maxmua])
   legend
end

function make_figure(PD,mcml_data_d,lambda)
%make figure with all skin Layers and the PD DO plot
    figure
    hold on
    Fz_size = size(mcml_data_d.Fz,1);
    PD = PD.*mcml_data_d.dz;
    plot(lambda,PD,'DisplayName','PD');
    depth = 0;
    for i = 1:size(mcml_data_d.d)
        depth = depth + mcml_data_d.d(i); 
        name = sprintf('Layer %d', i);
        plot([lambda(1) lambda(end)],[depth depth],'DisplayName',name);
    end
    xlabel('wavelength [nm]')
    ylabel('PD')
    axis([lambda(1) lambda(end) 0 Fz_size*mcml_data_d.dz])
    set(gca, 'YDir','reverse');
    if(length(mcml_data_d.d) == 6)
        title("Penetration Depth and Depth Origin for normal skin")
    else
        title("Penetration Depth and Depth Origin for compressed skin")
    end
    hold off
    legend
end