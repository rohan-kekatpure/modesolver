% This script generates the figure for the paper. It is assumed that the
% microsphere modes are already solved for with the following parameters:
% 
%       radius = 2.5 microns
%       delta_r = delta_z = 25 nm
%       dpml_r = dpml_z = 1 um
%       lambda_TM ~ 1596.931 nm
%       lambda_TE ~ 1562.059 nm
%
% With these parameters, the TM mode is the first solution and the TE mode
% is the third. The purpose of this script is just to generate a pretty
% redblue color-plot. 

function fplot_microsphere_modes_paper_figure(IN,OUT,R_sph)

close all;

um = 1e-6;

fs_xylabel = 20;
fs_ticklabel = 13;

ncont=300; % Set ncont to >200 for publication quality plots
x0 = 0.135;
y0 = 0.07;
right_clearence = 0.05;
top_clearence = 0.05;
hspacing = 0.05;
vspacing = 0.07;
W = (1 - x0 - hspacing - right_clearence)/2;
H = (1 - y0 - vspacing - top_clearence)/2;

r_i = IN.geom.r_i;
z_i = IN.geom.z_i;

nr = length(r_i);
nz = length(z_i);

Nzr = nr * nz;

r_min = IN.geom.r_min;
r_max = IN.geom.r_max;
z_min = IN.geom.z_min;
z_max = IN.geom.z_max;

% For some reason, the eigenvectors found by the present code are negative
% of the eigenvectors found by the previous code. While there is nothing
% wrong mathematically, I have to artificially make the values negative in
% order to make the EPS figures consistent with the PDF figures that went
% out with the review.

Exr_TM = -OUT(1).Er;
Ezr_TM = -OUT(1).Ez;
 
Exr_TE = -OUT(3).Er;
Ezr_TE = -OUT(3).Ez;

f1=figure(1);

%--------------------------------------------------------------------------

f11 = axes( 'pos', [x0  y0 + H + vspacing W H]);        
contourf(r_i/um,z_i/um,real(Exr_TE),ncont,'LineStyle','none');         
axis equal;         
axis([r_min r_max z_min z_max]/um);                
draw_geom_sphere(R_sph/um);        
xlabel('$r\,(\mu\mathrm{m})$ ','Fontsize',fs_xylabel,'interpreter','latex');
ylabel('$z\,(\mu\mathrm{m})$','Fontsize',fs_xylabel,'interpreter','latex');
set(f11, 'fontsize',fs_ticklabel);
colormap redblue;
set(gca,'clim',[-0.004 +0.004]);

text(5.3,3.2,'$E_r^{\mathsf{TE}}$','fontsize',25,'interpreter','latex');
text(6, -3.5,'(a)','fontsize',20,'fontname','Arial');

%--------------------------------------------------------------------------

%f12 = axes( 'pos',[x0 + W + hspacing y0 + H + vspacing W H]);
f12 = axes( 'pos',[x0 y0 W H]);
contourf(r_i/um,z_i/um,real(Ezr_TE),ncont,'LineStyle','none');         
draw_geom_sphere(R_sph/um);  
xlabel('$r\,(\mu\mathrm{m})$ ','Fontsize',fs_xylabel,'interpreter','latex');
ylabel('$z\,(\mu\mathrm{m})$','Fontsize',fs_xylabel,'interpreter','latex');
set(f12, 'fontsize',fs_ticklabel);
axis equal;         
axis([r_min r_max z_min z_max]/um);
colormap redblue;
set(gca,'clim',[-0.0004 +0.0004]);

text(5.3,3.2,'$E_z^{\mathsf{TE}}$','fontsize',25,'interpreter','latex');
text(6, -3.5,'(b)','fontsize',20,'fontname','Arial');

%--------------------------------------------------------------------------

f13 = axes( 'pos',[x0 + W + hspacing y0 + H + vspacing W H]);
contourf(r_i/um,z_i/um,real(Exr_TM),ncont,'LineStyle','none');         
axis equal;         
axis([r_min r_max z_min z_max]/um);                
draw_geom_sphere(R_sph/um);        
xlabel('$r\,(\mu\mathrm{m})$ ','Fontsize',fs_xylabel,'interpreter','latex');
ylabel('$z\,(\mu\mathrm{m})$','Fontsize',fs_xylabel,'interpreter','latex');
set(f13, 'fontsize',fs_ticklabel);
colormap redblue;
set(gca,'clim',[-0.004 +0.004]);

text(5.2,3.2,'$E_r^{\mathsf{TM}}$','fontsize',25,'interpreter','latex');
text(6, -3.5,'(c)','fontsize',20,'fontname','Arial');

%--------------------------------------------------------------------------

f13 = axes( 'pos',[x0 + W + hspacing y0 W H]);
contourf(r_i/um,z_i/um,real(Ezr_TM),ncont,'LineStyle','none');         
draw_geom_sphere(R_sph/um);        
xlabel('$r\,(\mu\mathrm{m})$ ','Fontsize',fs_xylabel,'interpreter','latex');
ylabel('$z\,(\mu\mathrm{m})$','Fontsize',fs_xylabel,'interpreter','latex');
set(f13, 'fontsize',fs_ticklabel);
axis equal;         
axis([r_min r_max z_min z_max]/um);
colormap redblue;
set(gca,'clim',[-0.03 +0.03]);

text(5.2,3.2,'$E_z^{\mathsf{TM}}$','fontsize',25,'interpreter','latex');
text(6, -3.5,'(d)','fontsize',20,'fontname','Arial');

%--------------------------------------------------------------------------        

print_pdf_1 = 1;

if print_pdf_1 == 1
    set(f1,'PaperUnits','inches');
    set(f1, 'PaperPositionMode', 'manual');
    set(f1, 'PaperPosition', [-0.9 0 9 8]);
    set(f1,'PaperSize',[7.5 7.7]);   
    print(f1,'-dpdf','-r300','microsphere_mode_profiles.pdf');
    print(f1,'-deps','-r300','microsphere_mode_profiles.eps');
end

    

function draw_geom_sphere(r_sph)

rectangle('Curvature',[1,1],...
          'Position',[-r_sph,-r_sph,2*r_sph,2*r_sph],...
          'EdgeColor',[.1 .5 .5],...
          'linewidth',1);