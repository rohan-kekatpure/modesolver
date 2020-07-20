% This function takes the helmholtz operator matrix and calculates its
% eigenvalues using the 'eigs' function. The number of eigenvalues to find
% and the guess values for the eigen-wavelength is contained in the
% solparams field of the IN structure. 
%
% Going further, it also creates a 1 X nsol size output structure array
% where nsol is the number of eigenvalues found. Each element of the output
% structure array is a strcture containing following fields:
%
%           lambda0 : real eigen-wavelength of the Nth solution
%           Q       : Quality factor of the Nth solution
%           Ex      : The ny * nx Ex field profile for the Nth mode
%           Ez      : The ny * nx Ez field profile for the Nth mode
%
% The output structure can be used for plotting the fields etc. To plot the
% Er field of the 5th mode, for example, you can use:
%
%           surface(r_i/um,...
%                   z_i/um,...
%                   real(OUT(5).Er));
%                


function SOL = get_modes(IN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unpack the needed parameters from the input IN structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 nm = IN.const.nm;

[ny,nx] = size(IN.mater.epsr) ;

lambda0 = IN.solparams.lambda0;
neff_guess = IN.solparams.neff_guess;

nsol=IN.solparams.neigs ;

nyx = ny * nx ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k0 = 2 * pi/lambda0;

Xguess=speye(2*nyx);

evals_guess = (k0 * neff_guess)^2;

opts.tol = eps;

% m_solving = show_info('Solving...','Status');

[evecs,evals]=eigs(IN.solparams.M_helmholtz,...
                   Xguess,...
                   nsol,...
                   evals_guess,...
                   opts);

               
% close(m_solving);

% m_reshape = show_info('Creating Output Structure','Status');

beta_list = sqrt(diag(evals)); % Remember the eigenvalue solver returns beta^2

neff_list = beta_list/k0;

SOL = struct('neff',0,...     
              'Ex',cell(1,nsol),...
              'Ey',cell(1,nsol));


for sn = 1 : nsol



    Ex1D=evecs(1:nyx,sn); Ex_mat=reshape(Ex1D,nx,ny).';
    Ey1D=evecs(nyx+1:2*nyx,sn); Ey_mat=reshape(Ey1D,nx,ny)';

    SOL(sn).neff = neff_list(sn);
    SOL(sn).Ex = Ex_mat;
    SOL(sn).Ey = Ey_mat;
    

end

% close(m_reshape);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print the roots in a tabular form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

box_printf('List of complex Eigen-mode indices','*',50);

for ii = 1 : nsol
    fprintf('%d\t n_eff = %2.10f  %2.10f \n', ...
            ii,...
            real(SOL(ii).neff),...
            imag(SOL(ii).neff));
    fprintf('-----------------------------------------------------\n');    
end
