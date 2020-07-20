% This script generates the figure for the paper. It is assumed that the
% microdisk modes are already solved for with the following parameters:
% 
%       radius = 1.5 microns
%       m = 14        
%       delta_r = delta_z = 25 nm
%       dpml_r = dpml_z = 1 um
%
% With these parameters, the fundamental mode is the first eigenvalue
% returned 

function fplot_microdisk_modes_paper_figure(IN,OUT,solnum,r_disk,h_disk)

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
       
Exr = OUT(solnum).Er ;
Ezr = OUT(solnum).Ez ;
f1=figure(1);

%--------------------------------------------------------------------------

f1a = axes( 'pos', [x0  y0 + H + vspacing W H]);        
contourf(r_i/um,z_i/um,real(Exr),ncont,'LineStyle','none');         
axis equal;         
axis([r_min r_max z_min z_max]/um);                
draw_geom_rect(r_i,z_i,r_disk,h_disk)       
xlabel('$r\,(\mu\mathrm{m})$ ','Fontsize',fs_xylabel,'interpreter','latex');
ylabel('$z\,(\mu\mathrm{m})$','Fontsize',fs_xylabel,'interpreter','latex');
set(f1a, 'fontsize',fs_ticklabel);
colormap redblue;
set(f1a,'clim',[-0.06 +0.06]);

text(0.25,3.4,'$\mathsf{Re}(E_r)$','fontsize',15,'interpreter','latex','color','k');
text(0.25, -3.5,'(a)','fontsize',20,'fontname','Arial','color','k');

%--------------------------------------------------------------------------


f1b = axes( 'pos',[x0 y0 W H]);
contourf(r_i/um,z_i/um,log10(abs(real(Exr))),ncont,'LineStyle','none');         
draw_geom_rect(r_i,z_i,r_disk,h_disk) 
xlabel('$r\,(\mu\mathrm{m})$ ','Fontsize',fs_xylabel,'interpreter','latex');
ylabel('$z\,(\mu\mathrm{m})$','Fontsize',fs_xylabel,'interpreter','latex');
set(f1b, 'fontsize',fs_ticklabel);
axis equal;         
axis([r_min r_max z_min z_max]/um);
colormap redblue;
set(f1b,'clim',[-10 -2]);

text(0.25,3.4,'$\log\vert\mathsf{Re}(E_r)\vert$','fontsize',15,'interpreter','latex','color','k');
text(0.25, -3.5,'(b)','fontsize',20,'fontname','Arial','color','w');

%--------------------------------------------------------------------------

f1c = axes( 'pos',[x0 + W + hspacing y0 + H + vspacing W H]);
contourf(r_i/um,z_i/um,real(Ezr),ncont,'LineStyle','none');         
axis equal;         
axis([r_min r_max z_min z_max]/um);                
draw_geom_rect(r_i,z_i,r_disk,h_disk)       
xlabel('$r\,(\mu\mathrm{m})$ ','Fontsize',fs_xylabel,'interpreter','latex');
ylabel('$z\,(\mu\mathrm{m})$','Fontsize',fs_xylabel,'interpreter','latex');
set(f1c, 'fontsize',fs_ticklabel);
colormap redblue;
set(f1c,'clim',[-0.004 +0.0042]);

text(0.25,3.4,'$\mathsf{Re}(E_z)$','fontsize',15,'interpreter','latex','color','k');
text(0.25, -3.5,'(c)','fontsize',20,'fontname','Arial','color','k');

%--------------------------------------------------------------------------

f1d = axes( 'pos',[x0 + W + hspacing y0 W H]);
contourf(r_i/um,z_i/um,log10(abs(real(Ezr))),ncont,'LineStyle','none');         
draw_geom_rect(r_i,z_i,r_disk,h_disk)       
xlabel('$r\,(\mu\mathrm{m})$ ','Fontsize',fs_xylabel,'interpreter','latex');
ylabel('$z\,(\mu\mathrm{m})$','Fontsize',fs_xylabel,'interpreter','latex');
set(f1d, 'fontsize',fs_ticklabel);
axis equal;         
axis([r_min r_max z_min z_max]/um);
colormap redblue;
set(f1d,'clim',[-10 -2]);

text(0.25,3.4,'$\log\vert\mathsf{Re}(E_z)\vert$','fontsize',15,'interpreter','latex','color','k');
text(0.25, -3.5,'(d)','fontsize',20,'fontname','Arial','color','w');

%--------------------------------------------------------------------------        

print_pdf_1 = 1;

if print_pdf_1 == 1
    set(f1,'PaperUnits','inches');
    set(f1, 'PaperPositionMode', 'manual');
    set(f1, 'PaperPosition', [-0.9 0 9 8]);
    set(f1,'PaperSize',[7.5 7.7]);   
    print(f1,'-dpdf','-r300','microdisk_mode_profiles.pdf');
    print(f1,'-deps','-r300','microdisk_mode_profiles.eps');
end    

function draw_geom_rect(r_i,z_i,r_disk,h_disk)

um = 1e-6;

delta_z = z_i(2) - z_i(1);

delta_r = r_i(2) - r_i(1);

rectangle('Position',[r_i(1),...
                      -h_disk/2 - delta_z/2,...
                      r_disk - r_i(1) - 3*delta_r/2,...
                      h_disk]/um,...
          'LineWidth',1.5,...
          'EdgeColor',[0.1,0.5,0.5]);
      
      
      