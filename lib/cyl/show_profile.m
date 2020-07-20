% This function generates a two-dimensional plot of either the relative
% permittivity (opt = 'index') or the PML's (opt = 'pml'). This functions
% is called from the main code. Once the epsilon profile is ensured to be
% correct, a call to this function is usually commented out.

function show_profile(IN,opt)

um = IN.const.um;


figure(1); 

switch opt
    case 'pml'
        
        subplot(1,2,1); 
        surface(IN.geom.r_i/um,IN.geom.z_i/um,-imag(IN.pml.sz_i),'linestyle','none'); 
        colorbar;
        title('z-pml profile', 'interpreter','latex','fontsize',15);
        axis([IN.geom.r_min,...
              IN.geom.r_max,...
              IN.geom.z_min,...
              IN.geom.z_max]/um,'square');

        subplot(1,2,2); 
        surface(IN.geom.r_i/um,IN.geom.z_i/um,-imag(IN.pml.sr_i),'linestyle','none'); 
        colorbar;
        title('r-pml profile', 'interpreter','latex','fontsize',15);
        axis([IN.geom.r_min,...
              IN.geom.r_max,...
              IN.geom.z_min,...
              IN.geom.z_max]/um,'square');
        


        colormap bluewhite

    case 'index'
        
        surface(IN.geom.r_i/um,IN.geom.z_i/um,IN.mater.epsr,'LineStyle','-');
        axis([IN.geom.r_min,...
              IN.geom.r_max,...
              IN.geom.z_min,...
              IN.geom.z_max]/um,'square');
        
    otherwise
        
        fprintf('\n\nUNKNOWN OPTION\n\n');
     
end

return;