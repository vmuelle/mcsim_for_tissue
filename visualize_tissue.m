function visualize_tissue(d_s,n_c,p)
   tissue = Tissue();
   tl = tissue.makeTissueList(d_s,n_c,p);

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
   ylabel('\mu_a [cm^{-1}]')
   axis([lambda(1) lambda(end) 0 maxmua])
   legend
end

