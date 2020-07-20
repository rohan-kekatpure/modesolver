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
        surface(IN.geom.x_i/um,IN.geom.y_i/um,-imag(IN.pml.Sy_ii),'linestyle','none'); 
        colorbar;
        title('y-pml profile', 'interpreter','latex','fontsize',15);
        axis([IN.geom.x_min,...
              IN.geom.x_max,...
              IN.geom.y_min,...
              IN.geom.y_max]/um,'square');

        subplot(1,2,2); 
        surface(IN.geom.x_i/um,IN.geom.y_i/um,-imag(IN.pml.Sx_ii),'linestyle','none'); 
        colorbar;
        title('x-pml profile', 'interpreter','latex','fontsize',15);
        axis([IN.geom.x_min,...
              IN.geom.x_max,...
              IN.geom.y_min,...
              IN.geom.y_max]/um,'square');
        
        colormap(MakeColorMap([1 1 1], [0 0 1], 25))

    case 'index'
        
        surface(IN.geom.x_i/um,IN.geom.y_i/um,IN.mater.epsr,'LineStyle','-');
        axis([IN.geom.x_min,...
              IN.geom.x_max,...
              IN.geom.y_min,...
              IN.geom.y_max]/um,'square');
        
    otherwise
        
        fprintf('\n\nUNKNOWN OPTION\n\n');
     
end

return;