close all force;
clear all ; 
clc;

format long;

% Cartesian mode solver
addpath("../lib/cart");


tic;
I = sqrt(-1);
um = 1e-6;
nm = 1e-9;
Z0 = 377;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometrical and material parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps_core = 2.0^2 ;
eps_clad = 1.45^2;
fradx = 500 * nm ;
frady = 525 * nm ;
dpml_x = 1 * um;
dpml_y = 1 * um;
M_pml = 2; % quadratic PML

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda0 = 1550 * nm ; % Wavelength for mode solution 
neff_guess = 1.8 - 0i; % Guess effective index of the sought mode
nsol = 6; % number of eigenmodes to find 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the mesh and co-ordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_x = 50 * nm;
x_min = -3 * um;
x_max = 3 * um;
delta_y = delta_x;                        
y_min = -3 * um;
y_max = +3 * um;


[x_i, x_h, y_i, y_h] = get_xy(x_min,x_max,delta_x,...
                              y_min,y_max,delta_y);

nx = numel(x_i);
ny = numel(y_i);

box_printf(['NUMBER OF X GRID POINTS = ', num2str(nx),...
            '\n\n',...
            'NUMBER OF Y GRID POINTS = ', num2str(ny)],'*',50);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating the epsilon profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


epsr = zeros(ny,nx); % nz rows, nr columns

for i = 1 : ny    
    for j = 1 : nx        
        if    ( (x_i(j)/fradx)^2 + (y_i(i)/frady)^2 ) < 1                 
                epsr(i,j) = eps_core;                
        else                 
                epsr(i,j) = eps_clad;
        end
    end
end

[epsr_xh, epsr_yh, epsr_xyhh] = get_shift_epsr(epsr);


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
IN.solparams.lambda0 = lambda0;
IN.solparams.neff_guess = neff_guess ;
IN.solparams.neigs = nsol ;
IN.solparams.pmlopt = 'off'; 

% Legal options are 'off' and 'on'. For faster and stabler solutions keep
% it  'off' unless you are looking specifically for leaky modes. Even then,
% you should first solve the problem using  the 'off' option and note down
% index of the mode you are looking for (say nx). Then turn 'on' the PML
% and provide nx as the initial guess.

                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the GEOMETRY field of the input structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IN.geom.x_min = x_min ;
IN.geom.x_max = x_max ;
IN.geom.y_min = y_min ;
IN.geom.y_max = y_max ;
IN.geom.delta_x = delta_x ;
IN.geom.delta_y = delta_y ;
IN.geom.x_i = x_i ;
IN.geom.x_h = x_h;
IN.geom.dpml_x = dpml_x; 
IN.geom.y_i = y_i ;
IN.geom.y_h = y_h ;
IN.geom.dpml_y = dpml_y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the EPSILON field of the input structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IN.mater.epsr = epsr;
IN.mater.epsr_xh = epsr_xh;
IN.mater.epsr_yh = epsr_yh;
IN.mater.epsr_xyhh = epsr_xyhh;
                        
box_printf(['THE FIBER HORIZONTAL DIAMETER SPANS ',...
               num2str(nnz(epsr(:,floor(nx/2)) == eps_core)),...
               ' GRID POINTS\n\n',...
              'THE FIBER VERTICAL DIAMETER SPANS ',...
               num2str(nnz(epsr(floor(ny/2),:) == eps_core)),...
               ' GRID POINTS'],'*',50);

[Sx_ii, Sx_hh, Sx_hi, Sy_ii, Sy_hh, Sy_ih] = get_pmls(IN,M_pml) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the PML field of the input structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IN.pml.Sx_ii = Sx_ii ; 
IN.pml.Sx_hi = Sx_hi ;
IN.pml.Sx_hh = Sx_hh ; 

IN.pml.Sy_ii = Sy_ii ; 
IN.pml.Sy_ih = Sy_ih ;
IN.pml.Sy_hh = Sy_hh ; 


show_profile(IN,'index');
% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLUTION of the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IN.solparams.M_helmholtz = get_helmholtz_mat(IN);
nnz(IN.solparams.M_helmholtz);
OUT = get_modes(IN);

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all force;