% This function organizes the raw eigenvlaues and eigenvectors returned in
% the solution step. It creates an array 'SOL' of solutions. Each element
% of 'SOL' possesses the following fields:
%
%                   lambda0 : real eigen-wavelength of the Nth solution
%                   Q       : Quality factor of the Nth solution
%                   Ex      : The nz * nr Ex field profile for the Nth mode
%                   Ez      : The nz * nr Ez field profile for the Nth mode
%
%

function SOL = get_field_profiles(IN,lambda_list,evecs)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unpack the needed parameters from the input IN structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_i = IN.geom.r_i;
z_i = IN.geom.z_i;

nsol = IN.solparams.neigs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now process the parameters as needed for this script 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nr = length(r_i);
nz = length(z_i);
Nzr = nr * nz;

Q_list = -real(lambda_list)./ (2 * imag(lambda_list));

%SOL = zeros(size(lambda_list));

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