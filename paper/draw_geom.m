% This is a customized function to draw microring geometry. It is called by
% the 'plotPDF_microring_modes_new.m'script

function draw_geom(IN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Unpack the 'IN' structure to extract ring parameters
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

um = IN.const.um;

r_min = IN.geom.r_min ; 
r_max = IN.geom.r_max ; 

R_ring = IN.geom.ring_params(1);
WT = IN.geom.ring_params(2);
WB = IN.geom.ring_params(3);
ZT = IN.geom.ring_params(4);
ZB = IN.geom.ring_params(5);

rc = [R_ring-WT/2 R_ring-WB/2 R_ring+WB/2 R_ring+WT/2 R_ring-WT/2]/um;
zc = [ZT ZB ZB ZT ZT]/um;

line(rc,zc,...   
      'Color',[.1 .5 .5],...
      'linewidth',1.);
  
line([r_min r_max]/um,[ZB ZB]/um,...
      'Color',[.1 .5 .5],...
      'linewidth',1.);

