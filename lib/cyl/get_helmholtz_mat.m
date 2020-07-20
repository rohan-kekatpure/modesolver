% get_helmholtz(input_structure) generates the vector helmholtz operator
% whose eigenvalues are the eigenfrequencies we seek. This functions is
% called from the main program after epsilon_profile, and the PMLs are
% created. 

function M_hemlholtz = get_helmholtz_mat(IN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unpack the needed parameters from the input IN structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mphi = IN.solparams.mphi;


delta_r = IN.geom.delta_r;

delta_z = IN.geom.delta_z;




epsr = IN.mater.epsr;

epsr_rh = IN.mater.epsr_rh;

epsr_zm = IN.mater.epsr_zm;

epsr_zp = IN.mater.epsr_zp;




rs_i = IN.pml.rs_i ; 

rs_h = IN.pml.rs_h ; 

sr_i = IN.pml.sr_i ;

sr_h = IN.pml.sr_h ;

sz_i = IN.pml.sz_i ;

sz_h = IN.pml.sz_h ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dr2 = delta_r^2;

dz2 = delta_z^2;

drdz = delta_r * delta_z;

[nz , nr] = size(epsr);

Nzr = nz*nr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of eigenvalue matrices : Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mrr1 = spalloc(Nzr,Nzr,3*Nzr);

Mrr2 = spalloc(Nzr,Nzr,3*Nzr);

Mrz1 = spalloc(Nzr,Nzr,3*Nzr);

Mrz2 = spalloc(Nzr,Nzr,3*Nzr);

Mzr1 = spalloc(Nzr,Nzr,3*Nzr);

Mzr2 = spalloc(Nzr,Nzr,3*Nzr);

Mzz2 = spalloc(Nzr,Nzr,3*Nzr);

Mzz1 = spalloc(Nzr,Nzr,3*Nzr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Mrr1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


wbmh=waitbar(0,'Constructing Helmholtz Operator...');

clear m n;

for n = 1 : nz
    
    for m = 1 : nr        
        
        % COEFFICIENT OF Er(n-1,m)        
            
        if n > 1
            
            elem_mrr1_1 = (...
                              1/(sz_h(n-1,m) * sz_h(n-1,m))...
                           )/dz2;
            
            Mrr1( (n-1)*nr + m, (n-2)*nr + m) = elem_mrr1_1/epsr_rh(n,m);
            
        end
        
        % COEFFICIENT OF Er(n,m)
        
        if n > 1
            elem_mrr1_2 = -(... 
                                1/sz_h(n,m)^2 + 1/(sz_h(n-1,m) * sz_h(n-1,m))...
                           )/dz2...                            
                           - 0*mphi^2/rs_h(n,m)^2;
            
            Mrr1( (n-1)*nr + m, (n-1)*nr + m) = elem_mrr1_2/epsr_rh(n,m);
            
        end
        
        % COEFFICIENT OF Er(n+1,m)
        
        if n < nz
            
            elem_mrr1_3 = (...
                            1/sz_h(n,m)^2 ...
                           )/dz2;
            
            Mrr1( (n-1)*nr + m, n*nr + m) = elem_mrr1_3/epsr_rh(n,m);
            
        end
        
    
    end  
     
    
end

waitbar(1/8,wbmh);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Mrr2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear m n; 


for n = 1 : nz
    
    for m = 1 : nr
        
        % COEFFICIENT OF Er(n,m-1)
        
        if m > 1
            
            elem_mrr2_1 = (...
                                1/sr_i(n,m)^2 ...
                                * rs_h(n,m-1)/rs_h(n,m)...
                                * epsr_rh(n,m-1)/epsr(n,m)...
                            )/dr2;
                        
            Mrr2((n-1)*nr+m,(n-1)*nr+m-1) = elem_mrr2_1/epsr_rh(n,m);
            
        end

        % COFFICIENT OF Er(n,m)
        
        if m < nr
            
            elem_mrr2_2 = -(...
                              1/( sr_i(n,m) * sr_i(n,m+1) )* epsr_rh(n,m)/epsr(n,m)...
                              +...
                              1/sr_i(n,m)^2 * rs_i(n,m+1)/rs_i(n,m) * epsr_rh(n,m)/epsr(n,m+1)...
                           )/dr2 ...
                           - mphi^2/rs_h(n,m)^2;

            Mrr2((n-1)*nr+m,(n-1)*nr+m) = elem_mrr2_2/epsr_rh(n,m);
        end

        %COEFFICIENT OF Er(n,m+1)
        
        if m < nr
            
            elem_mrr2_3 = (...
                            1/( sr_i(n,m) * sr_i(n,m+1) )...
                            * rs_i(n,m+1) * rs_h(n,m+1)/(rs_i(n,m) * rs_h(n,m))...
                            * epsr_rh(n,m+1)/epsr(n,m+1)...
                         )/dr2 ;

            Mrr2((n-1)*nr+m,(n-1)*nr+m+1) = elem_mrr2_3/epsr_rh(n,m);

        end
        
    end
    
    
end

waitbar(2/8,wbmh);

Mrr = -Mrr1 - Mrr2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Mrz1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear m n;



for n = 1 : nz
    
    for m = 1 : nr
        
        % COEFFICIENT OF Ez(n,m)
        
        elem_mrz1_1 = -(...
                           1/( sr_h(n,m) * sz_i(n,m) )...
                       )/drdz ; 
        
        Mrz1((n-1)*nr+m,(n-1)*nr+m) = elem_mrz1_1 / epsr_rh(n,m);
        
        if m < nr
            
            elem_mrz1_2 = (...
                                1/( sr_h(n,m) * sz_i(n,m) )...
                          )/drdz;
            
            Mrz1((n-1)*nr+m,(n-1)*nr+m+1) = elem_mrz1_2 / epsr_rh(n,m);
            
        end
        
        if n >= 2
            
            elem_mrz1_3 = (...
                                1/( sr_h(n,m) * sz_i(n,m) )...
                          )/drdz;
            
            Mrz1((n-1)*nr+m,(n-2)*nr+m) = elem_mrz1_3 / epsr_rh(n,m);
            
            elem_mrz1_4 = -(...
                                1/( sr_h(n,m) * sz_i(n,m) )...
                           )/drdz;
            
            Mrz1((n-1)*nr+m,(n-2)*nr+m+1) = elem_mrz1_4 / epsr_rh(n,m);
            
        end
        
    end 
    
end

waitbar(3/8,wbmh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Mrz2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear m n;



for n = 1 : nz
    
    for m = 1 : nr
        
        % COEFFICIENT OF Ez(n,m)
        
        elem_mrz2_1 = (...
                          -rs_i(n,m)/rs_h(n,m) ...
                        * epsr_zp(n,m)/epsr(n,m)...
                      )...
                      * ( 1/(sz_i(n,m) * sr_h(n,m)) )/drdz;
        
        Mrz2((n-1)*nr+m,(n-1)*nr+m) = elem_mrz2_1 / epsr_rh(n,m); 
        
        
        %COEFICIENT OF Ez(n,m+1)        
        
        if m < nr
            
            elem_mrz2_2 = ( ...
                                rs_i(n,m+1)^2/(rs_h(n,m) * rs_i(n,m))...
                                * epsr_zp(n,m+1)/epsr(n,m+1)...
                            )...
                            * ( 1/(sz_i(n,m) * sr_h(n,m)) )/drdz;
            
            Mrz2((n-1)*nr+m,(n-1)*nr+m+1) = elem_mrz2_2 / epsr_rh(n,m);
            
        end
        
        
        %COEFFICIENT OF Ez(n-1,m)
        
        if n >= 2
            
            elem_mrz2_3 = (...
                                    rs_i(n,m)/rs_h(n,m)...
                                   * epsr_zp(n-1,m)/epsr(n,m)...
                              )...
                              * ( 1/(sz_i(n,m) * sr_h(n,m)) )/drdz;
            
            Mrz2((n-1)*nr+m,(n-2)*nr+m) = elem_mrz2_3 / epsr_rh(n,m);
            
        end
        
        
        %COEFFICIENT OF Ez(n-1,m+1)
        
        if ( n >= 2 ) && (m < nr)
                
                elem_mrz2_4 = (...
                                    -rs_i(n,m+1)^2/(rs_h(n,m) * rs_i(n,m))...
                                  * epsr_zp(n-1,m+1)/epsr(n,m+1)...
                                )...
                                * ( 1/(sz_i(n,m) * sr_h(n,m)) )/drdz;
                
                Mrz2((n-1)*nr+m,(n-2)*nr+m+1) = elem_mrz2_4 / epsr_rh(n,m);
                
        end
        
    end
    
end

waitbar(4/8,wbmh);

Mrz = Mrz1 - Mrz2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Mzr1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear m n;



for n = 1 : nz
    
    for m = 1 : nr
        
        %COEFFICIENT OF Er(n,m)
        
        elem_mzr1_1 = (...
                        -1/(sz_h(n,m) * sr_i(n,m))...
                        )/drdz; 
        
        Mzr1( (n-1) * nr + m, (n-1) * nr + m) = elem_mzr1_1 / epsr_zp(n,m);
        
        
        % COEFFICIENT OF Er(n,m-1) 
        
        if m > 1
            
            elem_mzr1_2 = (...                               
                               rs_i(n,m-1)/rs_i(n,m)...
                            )...
                            * ( 1/(sz_h(n,m) * sr_i(n,m)) )/drdz;
            
            Mzr1( (n-1) * nr + m, (n-1) * nr + m - 1) = elem_mzr1_2 / epsr_zp(n,m);
            
        end
      
        % COEFFICIENT OF Er(n+1,m) 
        
        if n < nz
            
            elem_mzr1_3 =  (...
                            1/( sz_h(n,m) * sr_i(n,m) )...
                            )/drdz;
            
            Mzr1( (n-1) * nr + m, n * nr + m ) = elem_mzr1_3 / epsr_zp(n,m);
            
        end
        
        %COEFFICIENT OF Er(n+1,m-1)
        
        if (n < nz) && (m > 1)
            
            elem_mzr1_4 = (...
                                -rs_i(n,m-1)/rs_i(n,m)...
                             )...
                             * ( 1/( sz_h(n,m) * sr_i(n,m) ) )/drdz;
                        
            Mzr1( (n-1) * nr + m, n * nr + m - 1 ) = elem_mzr1_4 / epsr_zp(n,m);
            
        end
            
    end

end


waitbar(5/8,wbmh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Mzr2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear m n;



for n = 1 : nz
    
    for m = 1 : nr
        
        % COEFFICIENT OF Er(n,m)
        
        elem_mzr2_1 = (...
                           - rs_h(n,m)/rs_i(n,m)...
                           * epsr_rh(n,m)/epsr(n,m)...
                         )...
                         * ( 1/( sz_h(n,m) * sr_i(n,m) ) )/drdz;
        
        Mzr2( (n-1) * nr + m, (n-1) * nr + m) =  elem_mzr2_1 / epsr_zp(n,m);
        
        % COEFFICIENT OF Er(n,m-1)
        
        if m > 1
            
            elem_mzr2_2 = (...
                                rs_h(n,m-1)/rs_i(n,m) ...                
                                * epsr_rh(n,m-1)/epsr(n,m)...
                             )...
                             * ( 1/( sz_h(n,m) * sr_i(n,m) ) )/drdz;
            
            Mzr2( (n-1) * nr + m, (n-1) * nr + m - 1) = elem_mzr2_2 / epsr_zp(n,m);
            
        end
        
        % COEFFICIENT OF Er(n+1,m)
        
        if n < nz
            
            elem_mzr2_3 = (...
                                    rs_h(n,m)/rs_i(n,m)...
                                    * epsr_rh(n+1,m)/epsr(n+1,m)...
                             )...
                             * ( 1/( sz_h(n,m) * sr_i(n,m) ) )/drdz;
            
            Mzr2( (n-1) * nr + m, n * nr + m ) = elem_mzr2_3 / epsr_zp(n,m);
            
        end
                
        % COEFFICIENT OF Er(n+1,m-1)
        
        if (n < nz) && (m > 1)
            
            elem_mzr2_4 = (...
                                  -rs_h(n,m-1)/rs_i(n,m) ...
                                  * epsr_rh(n+1,m-1)/epsr(n+1,m)...
                             )...
                             * ( 1/( sz_h(n,m) * sr_i(n,m) ) )/drdz;
            
            Mzr2( (n-1) * nr + m, n * nr + m - 1 ) = elem_mzr2_4 / epsr_zp(n,m);
            
        end
            
    end

end

waitbar(6/8,wbmh);

Mzr = Mzr1 - Mzr2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Mzz1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear m n;


for n = 1 : nz
    
    for m = 1 : nr
        
        % COEFFICIENT OF Ez(n,m)
    
        if n < nz
            
            elem_mzz1_1 = -( epsr_zp(n,m)/epsr(n,m) * 1/( sz_i(n,m) * sz_i(n+1,m) )...
                                   + epsr_zp(n,m)/epsr(n+1,m) * 1/sz_i(n,m)^2 ...
                            )/dz2 ...
                            - mphi^2/rs_i(n,m)^2;
            
            Mzz1( (n-1)*nr + m, (n-1)*nr + m) = elem_mzz1_1/epsr_zp(n,m);
                                            
        end
        
        % COEFFICIENT OF Ez(n+1,m)
        
        if n < nz
            
            elem_mzz1_2 = (...
                                 epsr_zp(n+1,m)/epsr(n+1,m) * 1/( sz_i(n,m) * sz_i(n+1,m) )...
                                )/dz2;
            
            Mzz1( (n-1)*nr + m, n*nr + m) = elem_mzz1_2 / epsr_zp(n,m);                                
            
        end
        
        % COEFFICIENT OF Ez(n-1,m)
        
        if n > 1
            
            elem_mzz1_3 = (...
                                 epsr_zp(n-1,m)/epsr(n,m) * 1/sz_i(n,m)^2 ...
                                )/dz2; 
            
            Mzz1( (n-1)*nr + m, (n-2)*nr + m) = elem_mzz1_3 / epsr_zp(n,m);
            
        end
            
    end

end

waitbar(7/8,wbmh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of Mzz2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%wbmzz2=waitbar(0,'Constructing Mzz2');

clear m n;



for n = 1 : nz
    
    for m = 1 : nr
        
        % COEFFICIENT OF Er(n,m)
    
        if m > 1            
            
            elem_mzz2_1 = -( ...
                                  1/( sr_h(n,m) * sr_h(n,m-1) )...  
                                  + rs_i(n,m-1)/rs_i(n,m) * 1/sr_h(n,m)^2 ...
                             )/dr2 - 0*mphi^2/rs_i(n,m)^2;
            
            Mzz2( (n-1)*nr + m, (n-1)*nr + m) = elem_mzz2_1 / epsr_zp(n,m);
            
        end
        
        % COEFFICIENT OF Er(n,m-1)
        
        if m > 1
            
            elem_mzz2_2 =  (...
                                   rs_i(n,m-1)/rs_i(n,m)...
                                   * 1/( sr_h(n,m) * sr_h(n,m-1) )...
                              )/dr2;
            
            Mzz2( (n-1)*nr + m, (n-1)*nr + m-1) = elem_mzz2_2 / epsr_zp(n,m);
            
        end
        
        % COEFFICIENT OF Er(n,m+1)
        
        if m < nr
            
            elem_mzz2_3 = (...
                            1/sr_h(n,m)^2 ...
                            )/dr2 ;
            
            Mzz2( (n-1)*nr + m, (n-1)*nr + m+1) = elem_mzz2_3 / epsr_zp(n,m);
            
        end
            
    end
    
end

waitbar(8/8,wbmh);

Mzz =  -Mzz1 -Mzz2 ;

M_hemlholtz = [Mrr Mrz ; 
    
             Mzr Mzz  ];
 
close(wbmh); 