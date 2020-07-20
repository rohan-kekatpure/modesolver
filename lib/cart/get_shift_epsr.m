% This function takes the co-ordinates and the permittivity at integer
% grid-points and outputs the the following quantities:
%
%         r_h : radial co-ordinates of half-integer grid points
% 
%         z_h : z co-ordinate of half-integer grid points
% 
%         epsr_rh : permittivity shifted by half grid point in r
% 
%         epsr_zp : permittivity shifted by plus half grid point in z
% 
%         epsr_zm : permittivity shifted by minus half grid point in z
%
% These quantities are subsequently added as fields to the input structure
% by the main code. 



function [epsr_xh, epsr_yh, epsr_xyhh] = get_shift_epsr(epsr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the needed parameters from the input IN structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ny,nx] = size(epsr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now process the parameters and generate the outputs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


epsr_xh = 0.5 * (epsr + ...
                      [epsr(:,2:nx) epsr(:,nx)]);


epsr_yh = 0.5 * ( epsr + ...
                            [epsr(2:ny,:) ; epsr(ny,:)]);
                        
epsr_xyhh = 0.5 * ( epsr_xh + epsr_yh);
