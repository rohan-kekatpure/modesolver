% This script is copied from microsphere_main.m and is used to calculate
% modes of a silicon nitride microsphere clad with oxide. Just didnt want
% to disturb the working example presented in microsphere_main1.m so I am
% copying the file before messing around.   

close all;
clear all ; 
clc;

format long;

%tic;

I = sqrt(-1);

um = 1e-6;

nm = 1e-9;

Z0 = 377;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometrical and material parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps_sph = 2.5^2;

eps_amb = 1.45^2;

R_sph = 2.5 * um;

dpml_r = 1 * um;

dpml_z = 1 * um;

M_pml = 2; % quadratic PML

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_guess = (1484.93-0.0208235i)*nm; % Guess wavelength of the sought mode

mphi = 21;  % m-number of the sought mode

nsol = 20; % number of eigenvalues to find 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the mesh and co-ordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_r = 25 * nm;

r_min = delta_r;

r_max = 7 * um;

delta_z = delta_r ;                        
                        
z_min = -4 * um;

z_max = +4 * um;


[r_i, r_h, z_i, z_h] = get_rz(r_min,r_max,delta_r,...
                              z_min,z_max,delta_z);

nr = length(r_i);

nz = length(z_i);

box_printf(['NUMBER OF RADIAL GRID POINTS = ', num2str(nr),...
            '\n\n',...
            'NUMBER OF Z GRID POINTS = ', num2str(nz)],'*',50);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating the epsilon profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


epsr = zeros(nz,nr); % nz rows, nr columns

for i = 1 : nz
    
    for j = 1 : nr
        
        if (r_i(j)^2 + z_i(i)^2 <= R_sph^2)
                
                epsr(i,j) = eps_sph;
                
        else 
                
                epsr(i,j) = eps_amb;
        end
    end
end

[epsr_rh, epsr_zp,epsr_zm] = get_shift_epsr(epsr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% !!! DO NOT ALTER THE FOLLOWING CODE. 
%
% !!! ALL USER-DEFINED INPUT SHOULD BE DEFINED AT THIS TIME 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create the CONSTANTS and some of the SOLUTION parameters fields for the
% input structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


IN.const.I = I ;
IN.const.um = um ;
IN.const.nm = nm ;
IN.const.Z0 = Z0 ;
IN.solparams.lambda_guess = lambda_guess ;
IN.solparams.mphi = mphi ;
IN.solparams.neigs = nsol ;

                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the GEOMETRY field of the input structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IN.geom.r_min = r_min ;
IN.geom.r_max = r_max ;
IN.geom.z_min = z_min ;
IN.geom.z_max = z_max ;
IN.geom.delta_r = delta_r ;
IN.geom.delta_z = delta_z ;
IN.geom.r_i = r_i ;
IN.geom.r_h = r_h;
IN.geom.dpml_r = dpml_r; 
IN.geom.z_i = z_i ;
IN.geom.z_h = z_h ;
IN.geom.dpml_z = dpml_z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the EPSILON field of the input structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IN.mater.epsr = epsr;
IN.mater.epsr_rh = epsr_rh;
IN.mater.epsr_zm = epsr_zm;
IN.mater.epsr_zp = epsr_zp;
                        
box_printf(['THE MICROSPHERE LENGTH SPANS ',...
               num2str(nnz(epsr(:,1) == eps_sph)),...
               ' GRID POINTS\n\n',...
              'THE MICROSPHERE RADIUS SPANS ',...
               num2str(nnz(epsr(ceil(end/2),:) == eps_sph)),...
               ' GRID POINTS'],'*',50);

%  f1=surface(r_i/um,z_i/um,epsr_rh,'LineStyle','-');

[sr_i,sr_h,rs_i,rs_h,sz_i,sz_h] = get_pmls(IN,M_pml);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the PML field of the input structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IN.pml.sr_i = sr_i ;
IN.pml.sr_h = sr_h ;
IN.pml.rs_i = rs_i ;
IN.pml.rs_h = rs_h ;
IN.pml.sz_i = sz_i ;
IN.pml.sz_h = sz_h ;


%show_profile(IN,'index');
%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLUTION of the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IN.solparams.M_helmholtz = get_helmholtz_mat(IN);

tic;

OUT = get_modes(IN);

toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%