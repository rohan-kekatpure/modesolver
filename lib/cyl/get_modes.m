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
%           Ex      : The nz * nr Ex field profile for the Nth mode
%           Ez      : The nz * nr Ez field profile for the Nth mode
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

[nz,nr] = size(IN.mater.epsr) ;

lambda_guess = IN.solparams.lambda_guess ;

nsol=IN.solparams.neigs ;

Nzr = nz * nr ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xguess=speye(2*Nzr);

evals_guess = (2*pi./lambda_guess)^2;

m_solving = show_info('Solving...','Status');

[evecs,evals]=eigs(IN.solparams.M_helmholtz,...
                   Xguess,...
                   nsol,...
                   evals_guess);

               
close(m_solving);

m_reshape = show_info('Creating Output Structure','Status');

k0_list = sqrt(diag(evals)); % Remember the eigenvalue solver returns k0^2. 

lambda_list = 2*pi./k0_list;

Q_list = -real(lambda_list)./ (2 * imag(lambda_list));

SOL = struct('lambda0',0,...
              'Q',0,...
              'Er',cell(1,nsol),...
              'Ez',cell(1,nsol));


for sn = 1 : nsol



    Er1D=evecs(1:Nzr,sn); Er_mat=reshape(Er1D,nr,nz)';
    Ez1D=evecs(Nzr+1:2*Nzr,sn); Ez_mat=reshape(Ez1D,nr,nz)';

    SOL(sn).lambda0 = real(lambda_list(sn));
    SOL(sn).Q = Q_list(sn);
    SOL(sn).Er = Er_mat;
    SOL(sn).Ez = Ez_mat;
    

end

close(m_reshape);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print the roots in a tabular form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

box_printf('List of complex Eigen-wavelengths in nm','*',50);

for ii = 1 : nsol
    fprintf('%d\t lambda0 = %2.10f nm\t |\t Q = %2.10e \n', ...
            ii,...
            SOL(ii).lambda0/nm,...
            SOL(ii).Q);
    fprintf('-----------------------------------------------------');    
    fprintf('-----------------------------------------------------\n');    
end
