% This function generates PMLs on three sides of the epsilon profile. This
% cunction is called from the main function after the epsilon profile is
% generated and _before_ generating the vector helmholtz operator (by
% calling get_helmholtz_mat). 


function [sr_i,sr_h,rs_i,rs_h,sz_i,sz_h] = get_pmls(IN,M_pml)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unpack the needed parameters from the input IN structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = IN.const.I ;
um = IN.const.um ;
nm = IN.const.nm ;
Z0 = IN.const.Z0 ;
lambda_guess = IN.solparams.lambda_guess ;

r_min = IN.geom.r_min;
r_max = IN.geom.r_max;
z_min = IN.geom.z_min;
z_max = IN.geom.z_max;
delta_r = IN.geom.delta_r;
delta_z = IN.geom.delta_z;
r_i = IN.geom.r_i;
z_i = IN.geom.z_i;
r_h = IN.geom.r_h;
z_h = IN.geom.z_h;

dpml_r = IN.geom.dpml_r;
dpml_z = IN.geom.dpml_z;

epsr = IN.mater.epsr;
epsr_rh = IN.mater.epsr_rh;
epsr_zm = IN.mater.epsr_zm;
epsr_zp = IN.mater.epsr_zp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now process the parameters as needed for this script 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k0 = 2 * pi/lambda_guess;

[nz , nr] = size(epsr);

r_pml = r_max - dpml_r ;

z_pml_top = z_max - dpml_z ;

z_pml_bottom = z_min + dpml_z ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of PMLs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sr_i = ones(nz,nr);
sr_h = ones(nz,nr);
sz_i = ones(nz,nr);
sz_h = ones(nz,nr);

% We next define the radial stretching factors at integral (rsf_i) and
% half-integral (rsf_h) grid points. These stretching are unique to
% cylindrical PMLs. In the radial PML region, we replace any appearence of
% the radial co-ordinate 'r' with '(stretching_factor)*r' 

rsf_i = ones(nz,nr);
rsf_h = ones(nz,nr);

clear i j;
eta = 1; % put eta = 0 to turn off PMLs

sigma_max_r = (M_pml+1) / (150 * pi * delta_r);

sigma_max_z = (M_pml+1) / (150 * pi * delta_z);

for m = 1 : 1 : nr
    
    % CONSTRUCTION OF sr_i
    
    if r_i(m) >= r_pml
        
        sr_i(:,m) = 1 - ...
                  eta * I * Z0 * sigma_max_r./(k0 * epsr(:,m) ) * ...
                  ( (r_i(m) - r_pml) /dpml_r ) ^ M_pml ; 
        
        rsf_i(:,m) = 1 - eta * I * Z0 * sigma_max_r./( (M_pml+1) * k0 * epsr(:,m) ) * ...
                    ( (r_i(m) - r_pml) /dpml_r ) ^ (M_pml + 1);
              
    end
    
    % CONSTRUCTION OF sr_h
    
    if r_h(m) >= r_pml
        
        sr_h(:,m) = 1 - ...
                  eta * I * Z0 * sigma_max_r./(k0 * epsr_rh(:,m) ) * ...
                  ( (r_h(m) - r_pml) /dpml_r ) ^ M_pml ; 
        
        rsf_h(:,m) = 1 - eta * I * Z0 * sigma_max_r./( (M_pml+1) * k0 * epsr_rh(:,m) ) * ...
                    ( (r_h(m) - r_pml) /dpml_r ) ^ (M_pml + 1);
              
    end
    
end

for n = 1 : 1 : nz
    
    % Construction of sz_i
    
    if z_i(n) >= z_pml_top
        
        sz_i(n,:) = 1 - ...
                     eta * I * Z0 * sigma_max_z./(k0 * epsr(n,:) * dpml_z^M_pml) * ...
                  (z_i(n) - z_pml_top) ^ M_pml;          
        
    elseif z_i(n) <= z_pml_bottom
        
        sz_i(n,:) = 1 - ...
                     eta * I * Z0 * sigma_max_z./(k0 * epsr(n,:) * dpml_z^M_pml) * ...
                  (z_pml_bottom - z_i(n)) ^ M_pml;                 
        
    end
    
    % Construction of sz_h
    
    if z_h(n) >= z_pml_top
        
        sz_h(n,:) = 1 - ...
                     eta * I * Z0 * sigma_max_z./(k0 * epsr_zp(n,:) * dpml_z^M_pml) * ...
                  (z_h(n) - z_pml_top) ^ M_pml;          
        
    elseif z_h(n) <= z_pml_bottom
        
        sz_h(n,:) = 1 - ...
                     eta * I * Z0 * sigma_max_z./(k0 * epsr_zp(n,:) * dpml_z^M_pml) * ...
                  (z_pml_bottom - z_h(n)) ^ M_pml;                 
        
    end

    
end

rs_i = repmat(r_i,nz,1).*rsf_i;
rs_h = repmat(r_h,nz,1).*rsf_h;
