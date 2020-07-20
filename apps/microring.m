close all;
clear all ; 
clc;

format long;

path('/Users/rdkekat/Research/Codes/Matlab/ModeSolver_2D/FDFDcode/Cylindrical/FromFirstPrinciples/Modules',path);

I = sqrt(-1);

um = 1e-6;

nm = 1e-9;

Z0 = 377;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Geometrical and material parameters
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps_ring = 3.5^2; % permittivity of ring material

eps_air = 1.00^2; % permittivity of air

eps_substr = 1.45^2; % permittivity of the substrate

R_ring = 2.5 * um; % Average radius of the microring (at the midpoint)

h_ring = 250 * nm; % Thickness of the nitride layer defines ring thickness

WT = 300 * nm; % Top width of the trapezium

WB = 500 * nm; % Bottom width of the trapezium

ZB = 0; % Substrate sits on  z = 0 plane

ZT = ZB + h_ring; 

tan_theta = 0.5 * (WB - WT) / (ZT - ZB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Solution parameters
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lambda_guess = 1550 * nm; % Guess wavelength of the sought mode

lambda_guess = (1597.359 - I * 1.94e-11) * nm; % Guess wavelength of the sought mode

mphi = 24;  % m-number of the sought mode

nsol = 20; % number of eigenvalues to find 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% Defining the mesh and co-ordinates
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_r = 25 * nm;

r_min = delta_r;

r_max = 7 * um;

delta_z = delta_r ;                        
                        
z_min = -4 * um;

z_max = +4 * um;

dpml_r = 1 * um;

dpml_z = 1 * um;

M_pml = 2; % quadratic PML

[r_i, r_h, z_i, z_h] = get_rz(r_min,r_max,delta_r,...
                              z_min,z_max,delta_z);

nr = length(r_i);

nz = length(z_i);

box_printf(['NUMBER OF RADIAL GRID POINTS = ', num2str(nr),...
            '\n\n',...
            'NUMBER OF Z GRID POINTS = ', num2str(nz)],'*',50);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Creating the epsilon profile
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsr = ones(nz,nr); % nz rows, nr columns

for i = 1 : nz % iteration over z
    
    for j = 1 : nr % iteration over r
        
        if z_i(i) <= 0
            
            epsr(i,j) = eps_substr;
            
        end
        
        if (z_i(i) >= ZB) && (z_i(i) <=ZT)
            
            dw = (ZT - z_i(i)) * tan_theta;

            rR = R_ring + WT/2 + dw;  

            rL = R_ring - WT/2 - dw ; 
            
            if  (r_i(j) > rL) && (r_i(j) <= rR+eps)        
                    
                epsr(i,j) = eps_ring;
                
            end
                
        end
                
    end
    
end

[epsr_rh, epsr_zp,epsr_zm] = get_shift_epsr(epsr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% !!! !!! !!! !!! !!! !!! !!!  !!! STOP !!! !!! !!! !!! !!! !!! !!! !!! !!!
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
%
% Create the GEOMETRY field of the input structure
%
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

% "ring_params" field isnt needed by  any of the solver modules but is
% needed b needed by plotter function. We create this extra field so we
% dont have to pass the individual ring dimensions to the plotter function

IN.geom.ring_params = [R_ring WT WB ZT,ZB];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create the EPSILON field of the input structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IN.mater.epsr = epsr;
IN.mater.epsr_rh = epsr_rh;
IN.mater.epsr_zm = epsr_zm;
IN.mater.epsr_zp = epsr_zp;
                        
%  f1=surface(r_i/um,z_i/um,epsr_rh,'LineStyle','-');

[sr_i,sr_h,rs_i,rs_h,sz_i,sz_h] = get_pmls(IN,M_pml);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create the PML field of the input structure
%
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
%
% SOLUTION of the problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IN.solparams.M_helmholtz = get_helmholtz_mat(IN);

OUT = get_modes(IN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%