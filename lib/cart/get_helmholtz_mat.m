% get_helmholtz(input_structure) generates the vector helmholtz operator
% whose eigenvalues are the eigenfrequencies we seek. This functions is
% called from the main program after epsilon_profile, and the PMLs are
% created. 

function M_hemlholtz = get_helmholtz_mat(IN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unpack the needed parameters from the input IN structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda0 = IN.solparams.lambda0;


delta_x = IN.geom.delta_x;

delta_y = IN.geom.delta_y;




epsr = IN.mater.epsr;

epsr_xh = IN.mater.epsr_xh;

epsr_yh = IN.mater.epsr_yh;

epsr_xyhh = IN.mater.epsr_xyhh;


Sx_ii = IN.pml.Sx_ii ;
Sx_hi = IN.pml.Sx_hi ;
Sx_hh = IN.pml.Sx_hh ; 

Sy_ii = IN.pml.Sy_ii ;
Sy_ih = IN.pml.Sy_ih ;
Sy_hh = IN.pml.Sy_hh ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx2 = delta_x^2;

dy2 = delta_y^2;

dxdy = delta_x * delta_y;

[ny , nx] = size(epsr);

Nyx = ny * nx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of eigenvalue matxices : Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mxx1 = spalloc(Nyx,Nyx,3*Nyx);

Mxx2 = spalloc(Nyx,Nyx,3*Nyx);

Mxy1 = spalloc(Nyx,Nyx,3*Nyx);

Mxy2 = spalloc(Nyx,Nyx,3*Nyx);

Myx1 = spalloc(Nyx,Nyx,3*Nyx);

Myx2 = spalloc(Nyx,Nyx,3*Nyx);

Myy2 = spalloc(Nyx,Nyx,3*Nyx);

Myy1 = spalloc(Nyx,Nyx,3*Nyx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Mxx1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k0 = 2 * pi/lambda0;

wbmh=waitbar(0,'Constructing Helmholtz Operator...');

clear m n;

for n = 1 : ny
    
    for m = 1 : nx        
        
        % COEFFICIENT OF Er(n-1,m)        
            
        if n > 1
            
            Mxx1( (n-1)*nx + m, (n-2)*nx + m) = -1/Sy_hh(n-1,m)^2 ;
            
        end
        
        % COEFFICIENT OF Er(n,m)
        
        if n > 1
            
            Mxx1( (n-1)*nx + m, (n-1)*nx + m) = ...
                                2/(Sy_hh(n-1,m) * Sy_hh(n,m)) ...
                                - k0^2 * epsr_xh(n,m) * dy2 ;
            
            % We age going to divide Mxx1 by dy2 in the end. So the dy2
            % multiplying "k0^2 * epsr_xh * dy2" will drop down. This way
            % of doing things speeds up the Helmholtz assembly by having
            % the 'for' loops do fewer operations.            
                           
            
        end
        
        % COEFFICIENT OF Er(n+1,m)
        
        if n < ny
            
            Mxx1( (n-1)*nx + m, n*nx + m) = -1/Sy_hh(n,m)^2 ;
            
            
        end       
        
        
    end  
     
    
end

Mxx1 = Mxx1 / dy2 ;  % division by dy2 as promised

waitbar(1/8,wbmh);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Mxx2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear m n; 


for n = 1 : ny
    
    for m = 1 : nx
        
        % COEFFICIENT OF Er(n,m-1)
        
        if m > 1
            
            Mxx2((n-1)*nx+m,(n-1)*nx+m-1) = ...
                - epsr_xh(n,m-1)/epsr(n,m) * 1/Sx_ii(n,m)^2 ; 
                        
            
        end

        % COFFICIENT OF Er(n,m)
        
        if m < nx
            
            Mxx2((n-1)*nx+m,(n-1)*nx+m) = ...
                epsr_xh(n,m)/(Sx_ii(n,m+1) * Sx_ii(n,m)) ...
                * ( 1/epsr(n,m) + 1/epsr(n,m+1) );  
            
        end

        %COEFFICIENT OF Er(n,m+1)
        
        if m < nx
            
            Mxx2((n-1)*nx+m,(n-1)*nx + m+1) = ...
                - epsr_xh(n,m+1)/epsr(n,m+1) * 1/Sx_ii(n,m+1)^2 ; 
        end
        
    end    
    
end

Mxx2 = Mxx2 / dx2;

waitbar(2/8,wbmh);

Mxx = Mxx1 + Mxx2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Mxy1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear m n;



for n = 1 : ny
    
    for m = 1 : nx
        
        % COEFFICIENT OF Ez(n,m)
        
       Mxy1((n-1)*nx+m,(n-1)*nx+m) = - 1/(Sx_hh(n,m) * Sy_ii(n,m)); 
    
        
        % COEFFICIENT OF Ez(n,m+1)
        
        if m < nx
            
            Mxy1((n-1)*nx+m,(n-1)*nx+m+1) = 1/(Sx_hh(n,m) * Sy_ih(n,m));
            
        end
        
        % COEFFICIENT OF Ez(n-1,m)
                
        if n > 1
            
            
            Mxy1((n-1)*nx+m,(n-2)*nx+m) = 1/(Sx_hh(n-1,m) * Sy_ii(n,m));
            
        end    
            
        % COEFFICIENT OF Ez(n-1,m+1)
        
        if (n > 1) && (m < nx)
            
            Mxy1((n-1)*nx+m,(n-2)*nx+m+1) = - 1/(Sx_hh(n-1,m) * Sy_ih(n,m));
            
        end
        
    end 
    
end

Mxy1 = Mxy1 / dxdy;

waitbar(3/8,wbmh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Mxy2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear m n;



for n = 1 : ny
    
    for m = 1 : nx
        
        % COEFFICIENT OF Ez(n,m)
        
        Mxy2((n-1)*nx+m,(n-1)*nx+m) = ...
            epsr_yh(n,m)/epsr(n,m) * 1/(Sy_ii(n,m) * Sx_hh(n,m));
             
        
        
        %COEFICIENT OF Ez(n,m+1)        
        
        if m < nx
            
            Mxy2((n-1)*nx+m,(n-1)*nx+m+1) = ...
                - epsr_yh(n,m+1)/epsr(n,m+1) ...
                * 1/(Sy_ii(n,m+1) * Sx_hh(n,m)) ;
            
        end
        
        
        %COEFFICIENT OF Ez(n-1,m)
        
        if n >= 2
            
            Mxy2((n-1)*nx+m,(n-2)*nx+m) = ...
                - epsr_yh(n-1,m)/epsr(n,m) ...
                * 1/(Sy_ii(n,m) * Sx_hh(n-1,m)) ;

            
        end
        
        
        %COEFFICIENT OF Ez(n-1,m+1)
        
        if ( n >= 2 ) && (m < nx)
                
            Mxy2((n-1)*nx+m,(n-2)*nx+m+1) = ...
                 epsr_yh(n-1,m+1)/epsr(n,m+1) ...
                * 1/(Sy_ii(n,m+1) * Sx_hh(n-1,m)) ;
                    
                
        end
        
    end
    
end

Mxy2 = Mxy2 / dxdy;

Mxy = Mxy1 + Mxy2 ;

waitbar(4/8,wbmh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Myx1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear m n;



for n = 1 : ny
    
    for m = 1 : nx
        
        %COEFFICIENT OF Er(n,m)
        
        
        Myx1( (n-1) * nx + m, (n-1) * nx + m) = ...
            -1/(Sy_hh(n,m) * Sx_ii(n,m));
        
        
        % COEFFICIENT OF Er(n,m-1) 
        
        if m > 1
            
            Myx1( (n-1) * nx + m, (n-1) * nx + m - 1) = ...
                1/(Sy_hh(n,m-1) * Sx_ii(n,m));
            
        end
      
        % COEFFICIENT OF Er(n+1,m) 
        
        if n < ny
            
            Myx1( (n-1) * nx + m, n * nx + m ) = ...
                1/(Sy_hh(n,m) * Sx_hi(n,m));
            
        end
        
        %COEFFICIENT OF Er(n+1,m-1)
        
        if (n < ny) && (m > 1)
            
            Myx1( (n-1) * nx + m, n * nx + m - 1 ) = ...
                - 1/(Sy_hh(n,m-1) * Sx_hi(n,m));
            
        end
            
    end

end

Myx1 = Myx1 / dxdy;

waitbar(5/8,wbmh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Myx2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear m n;
 
for n = 1 : ny
    
    for m = 1 : nx
        
        % COEFFICIENT OF Er(n,m)
        
        Myx2( (n-1) * nx + m, (n-1) * nx + m) =...
            epsr_xh(n,m)/epsr(n,m) * 1/(Sx_ii(n,m) * Sy_hh(n,m)) ;
        
        % COEFFICIENT OF Er(n,m-1)
        
        if m > 1
            
            Myx2( (n-1) * nx + m, (n-1) * nx + m - 1) = ...
               - epsr_xh(n,m-1)/epsr(n,m) * 1/(Sx_ii(n,m) * Sy_hh(n,m-1)) ;
            
        end
        
        % COEFFICIENT OF Er(n+1,m)
        
        if n < ny
            
            Myx2( (n-1) * nx + m, n * nx + m ) = ...
                - epsr_xh(n+1,m)/epsr(n+1,m) * 1/(Sx_ii(n+1,m) * Sy_hh(n,m)) ;
            
        end
                
        % COEFFICIENT OF Er(n+1,m-1)
        
        if (n < ny) && (m > 1)
            
            Myx2( (n-1) * nx + m, n * nx + m - 1 ) = ...
                epsr_xh(n+1,m-1)/epsr(n+1,m) * 1/(Sx_ii(n+1,m) * Sy_hh(n,m-1)) ;
            
        end
            
    end

end

Myx2 = Myx2 / dxdy;
waitbar(6/8,wbmh);

Myx = Myx1 + Myx2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Myy1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear m n;


for n = 1 : ny
    
    for m = 1 : nx
        
        % COEFFICIENT OF Ez(n,m)
    
        if m > 1
            
            Myy1( (n-1)*nx + m, (n-1)*nx + m) = ...
                2/(Sx_hh(n,m-1) * Sx_hh(n,m));
                                            
        end
        
        % COEFFICIENT OF Ez(n,m+1)
        
        if m < nx
            
            Myy1( (n-1)*nx + m, (n-1)*nx + m+1) = - 1/Sx_hh(n,m)^2;                                
            
        end
        
        % COEFFICIENT OF Ez(n,m-1)
        
        if m > 1
            
            Myy1( (n-1)*nx + m, (n-1)*nx + m-1) = - 1/Sx_hh(n,m-1)^2;
            
        end
            
    end

end

Myy1 = Myy1 / dx2;

waitbar(7/8,wbmh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Myy2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%wbMyy2=waitbar(0,'Constructing Myy2');

clear m n;



for n = 1 : ny
    
    for m = 1 : nx
        
        % COEFFICIENT OF Er(n,m)
    
        if n < ny
            
            Myy2( (n-1)*nx + m, (n-1)*nx + m) = ...
                epsr_yh(n,m) / (Sy_ii(n,m) * Sy_ii(n+1,m))...
                * (1/epsr(n,m) + 1/epsr(n+1,m)) ...
                 - k0^2 * epsr_yh(n,m) * dy2 ;
            
        end
        
        % COEFFICIENT OF Er(n+1,m)
        
        if n < ny
            
 
            Myy2( (n-1)*nx + m, n*nx + m) = ...
                - epsr_yh(n+1,m)/epsr(n+1,m) * 1/Sy_ii(n+1,m)^2 ;
            
        end
        
        % COEFFICIENT OF Er(n-1,m)
        
        if n > 1
            
            Myy2( (n-1)*nx + m, (n-2)*nx + m) = ...
                - epsr_yh(n-1,m)/epsr(n,m) * 1/Sy_ii(n,m)^2 ;
            
        end
            
    end
    
end

Myy2 = Myy2 / dy2; 

waitbar(8/8,wbmh);

Myy =  Myy1 + Myy2 ;

M_hemlholtz = -[Mxx Mxy ; 
    
             Myx Myy  ];
 
close(wbmh); 