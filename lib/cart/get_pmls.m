% This function generates PMLs on all four sides of the epsilon profile.
% This function is called from the main function after the epsilon profile
% is generated and _before_ generating the vector helmholtz operator (by
% calling get_helmholtz_mat). 


function [Sx_ii, Sx_hh, Sx_hi, Sy_ii, Sy_hh, Sy_ih] = get_pmls(IN,M_pml)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unpack the needed parameters from the input IN structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = IN.const.I ;
um = IN.const.um ;
nm = IN.const.nm ;
Z0 = IN.const.Z0 ;
lambda0 = IN.solparams.lambda0;

x_min = IN.geom.x_min;
x_max = IN.geom.x_max;
y_min = IN.geom.y_min;
y_max = IN.geom.y_max;
delta_x = IN.geom.delta_x;
delta_y = IN.geom.delta_y;
x_i = IN.geom.x_i;
y_i = IN.geom.y_i;
x_h = IN.geom.x_h;
y_h = IN.geom.y_h;

dpml_x = IN.geom.dpml_x;
dpml_y = IN.geom.dpml_y;

epsr = IN.mater.epsr;
epsr_xh = IN.mater.epsr_xh;
epsr_yh = IN.mater.epsr_yh;
epsr_xyhh = IN.mater.epsr_xyhh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now process the parameters as needed for this script 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k0 = 2 * pi/lambda0;

[ny , nx] = size(epsr);

x_pml_right = x_max - dpml_x ;

x_pml_left = x_min + dpml_x ;

y_pml_top = y_max - dpml_y ;

y_pml_bottom = y_min + dpml_y ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of PMLs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sx_ii = ones(ny,nx);
Sx_hi = ones(ny,nx);
Sx_hh = ones(ny,nx);

Sy_ii = ones(ny,nx);
Sy_ih = ones(ny,nx);
Sy_hh = ones(ny,nx);

clear i j;
%eta = 1; % put eta = 0 to turn off PMLs <<ORIGINAL CODE>>

%<< MODIFICATION ON 8/29/11 >>

if strcmp(IN.solparams.pmlopt,'off') == 1
    eta = 0;
elseif strcmp(IN.solparams.pmlopt,'on') == 1
    eta = 1; 
else
    error('Wrong PML option. Legal options are ''off'' and ''on'' ');
end

sigma_max_x = (M_pml+1) / (150 * pi * delta_x);

sigma_max_y = (M_pml+1) / (150 * pi * delta_y);

for m = 1 : 1 : nx
    
    % CONSTRUCTION OF Sx_ii
    
    if x_i(m) >= x_pml_right
        
        Sx_ii(:,m) = 1 - ...
                  eta * I * Z0 * sigma_max_x./(k0 * epsr(:,m) ) * ...
                  ( (x_i(m) - x_pml_right) /dpml_x ) ^ M_pml ; 
      
    elseif x_i(m) <= x_pml_left
        
        Sx_ii(:,m) = 1 - ...
                  eta * I * Z0 * sigma_max_x./(k0 * epsr(:,m) ) * ...
                  ( ( x_pml_left - x_i(m)) /dpml_x ) ^ M_pml ; 
        
    end
              
    
    % CONSTRUCTION OF Sx_hi
    
    if x_i(m) >= x_pml_right
        
        Sx_hi(:,m) = 1 - ...
                  eta * I * Z0 * sigma_max_x./(k0 * epsr_yh(:,m) ) * ...
                  ( (x_i(m) - x_pml_right) /dpml_x ) ^ M_pml ;       
              
    elseif x_i(m) <= x_pml_left
        
        Sx_hi(:,m) = 1 - ...
                  eta * I * Z0 * sigma_max_x./(k0 * epsr_yh(:,m) ) * ...
                  ( ( x_pml_left - x_i(m) ) /dpml_x ) ^ M_pml ;       
    end
    
    % CONSTRUCTION OF Sx_hh
    
    if x_h(m) >= x_pml_right
        
        Sx_hh(:,m) = 1 - ...
                  eta * I * Z0 * sigma_max_x./(k0 * epsr_xyhh(:,m) ) * ...
                  ( (x_h(m) - x_pml_right) /dpml_x ) ^ M_pml ;       
              
    elseif x_h(m) <= x_pml_left
        
        Sx_hh(:,m) = 1 - ...
                  eta * I * Z0 * sigma_max_x./(k0 * epsr_xyhh(:,m) ) * ...
                  ( ( x_pml_left - x_h(m) ) /dpml_x ) ^ M_pml ;       
    end
 
    
end

% CONSTRUCTION OF Y-PMLS

for n = 1 : 1 : ny
    
    % Construction of Sy_ii
    
    if y_i(n) >= y_pml_top
        
        Sy_ii(n,:) = 1 - ...
                     eta * I * Z0 * sigma_max_y./(k0 * epsr(n,:) * dpml_y^M_pml) * ...
                  (y_i(n) - y_pml_top) ^ M_pml;          
        
    elseif y_i(n) <= y_pml_bottom
        
        Sy_ii(n,:) = 1 - ...
                     eta * I * Z0 * sigma_max_y./(k0 * epsr(n,:) * dpml_y^M_pml) * ...
                  (y_pml_bottom - y_i(n)) ^ M_pml;                 
        
    end

    % Construction of Sy_ih
    
    if y_i(n) >= y_pml_top
        
        Sy_ih(n,:) = 1 - ...
                     eta * I * Z0 * sigma_max_y./(k0 * epsr_xh(n,:) * dpml_y^M_pml) * ...
                  (y_i(n) - y_pml_top) ^ M_pml;          
        
    elseif y_i(n) <= y_pml_bottom
        
        Sy_ih(n,:) = 1 - ...
                     eta * I * Z0 * sigma_max_y./(k0 * epsr_xh(n,:) * dpml_y^M_pml) * ...
                  (y_pml_bottom - y_i(n)) ^ M_pml;                 
        
    end
    
    % Construction of Sy_hh
    
    if y_h(n) >= y_pml_top
        
        Sy_hh(n,:) = 1 - ...
                     eta * I * Z0 * sigma_max_y./(k0 * epsr_xyhh(n,:) * dpml_y^M_pml) * ...
                  (y_h(n) - y_pml_top) ^ M_pml;          
        
    elseif y_h(n) <= y_pml_bottom
        
        Sy_hh(n,:) = 1 - ...
                     eta * I * Z0 * sigma_max_y./(k0 * epsr_xyhh(n,:) * dpml_y^M_pml) * ...
                  (y_pml_bottom - y_h(n)) ^ M_pml;                 
        
    end

    
end
