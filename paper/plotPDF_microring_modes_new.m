% This function contains highly customized values for generating
% publishable PDF/EPS plots of microring modes. It is assumed that the main
% script 'microring_main.m' has been run and 20 eigen-frequencies have been
% output. You will have to mess around with the following command to set
% the correct value for the color limits:
%
%           set(gca,'clim',[xx,yy]
%
% Setting the correct values of xx and yy will allow you to get nice field
% plots. To get the figures centered on page with minimum of surrounding
% white space, play with the values of x0, y0, right_/top_clearence,
% and h/vspacing until you get what you want. 

function plotPDF_microring_modes_new(IN,OUT)

close all;

um = 1e-6;

fs_xylabel = 20;
fs_ticklabel = 13;


x0 = 0.135;
y0 = 0.07;
right_clearence = 0.05;
top_clearence = 0.05;
hspacing = 0.05;
vspacing = 0.07;
W = (1 - x0 - hspacing - right_clearence)/2;
H = (1 - y0 - vspacing - top_clearence)/2;

um = IN.const.um;

r_i = IN.geom.r_i;
z_i = IN.geom.z_i;

[nz nr] = size(IN.mater.epsr);

Nzr = nr * nz;

r_min = r_i(1);
r_max = r_i(nr);
z_min = z_i(1);
z_max = z_i(nz);

ncont=350; % Set ncont to >200 for publication quality plots
solnum = 2;

Er = OUT(solnum).Er ;
Ez = OUT(solnum).Ez ;

f1=figure(1);

%--------------------------------------------------------------------------

f11 = axes( 'pos', [x0  y0 + H + vspacing W H]);        
contourf(r_i/um,z_i/um,real(Er),ncont,'LineStyle','none');         
axis equal;         
axis([r_min r_max z_min z_max]/um);                
draw_geom(IN);        
xlabel('$r\,(\mu\mathrm{m})$ ','Fontsize',fs_xylabel,'interpreter','latex');
ylabel('$z\,(\mu\mathrm{m})$','Fontsize',fs_xylabel,'interpreter','latex');
set(f11, 'fontsize',fs_ticklabel);
colormap redblue;
set(gca,'clim',[-0.03 +0.03]);

text(5.3,3.2,'$E_r^{\mathsf{TE}}$','fontsize',25,'interpreter','latex');
text(6, -3.5,'(a)','fontsize',20,'fontname','Arial');

%--------------------------------------------------------------------------

%f12 = axes( 'pos',[x0 + W + hspacing y0 + H + vspacing W H]);
f12 = axes( 'pos',[x0 y0 W H]);
contourf(r_i/um,z_i/um,log10(abs(real(Er))),ncont,'LineStyle','none');         
draw_geom(IN);  
xlabel('$r\,(\mu\mathrm{m})$ ','Fontsize',fs_xylabel,'interpreter','latex');
ylabel('$z\,(\mu\mathrm{m})$','Fontsize',fs_xylabel,'interpreter','latex');
set(f12, 'fontsize',fs_ticklabel);
axis equal;         
axis([r_min r_max z_min z_max]/um);
colormap redblue;
set(gca,'clim',[-7 -1]);

text(2.8,3.2,'$\log_{10}(E_r^{\mathsf{TE}})$','fontsize',25,'interpreter','latex');
text(6, -3.5,'(b)','fontsize',20,'fontname','Arial');

%--------------------------------------------------------------------------

f13 = axes( 'pos',[x0 + W + hspacing y0 + H + vspacing W H]);
contourf(r_i/um,z_i/um,real(Ez),ncont,'LineStyle','none');         
axis equal;         
axis([r_min r_max z_min z_max]/um);                
draw_geom(IN);        
xlabel('$r\,(\mu\mathrm{m})$ ','Fontsize',fs_xylabel,'interpreter','latex');
ylabel('$z\,(\mu\mathrm{m})$','Fontsize',fs_xylabel,'interpreter','latex');
set(f13, 'fontsize',fs_ticklabel);
colormap redblue;
set(gca,'clim',[-0.03 0.03]);

text(5.3,3.2,'$E_z^{\mathsf{TE}}$','fontsize',25,'interpreter','latex');
text(6, -3.5,'(c)','fontsize',20,'fontname','Arial');

%--------------------------------------------------------------------------

f13 = axes( 'pos',[x0 + W + hspacing y0 W H]);
contourf(r_i/um,z_i/um,log10(abs(real(Ez))),ncont,'LineStyle','none');         
draw_geom(IN);        
xlabel('$r\,(\mu\mathrm{m})$ ','Fontsize',fs_xylabel,'interpreter','latex');
ylabel('$z\,(\mu\mathrm{m})$','Fontsize',fs_xylabel,'interpreter','latex');
set(f13, 'fontsize',fs_ticklabel);
axis equal;         
axis([r_min r_max z_min z_max]/um);
colormap redblue;
set(gca,'clim',[-7 -1]);

text(2.8,3.2,'$\log_{10}(E_z^{\mathsf{TE}})$','fontsize',25,'interpreter','latex');
text(6, -3.5,'(d)','fontsize',20,'fontname','Arial');

%--------------------------------------------------------------------------        

print_pdf_1 = 1;

if print_pdf_1 == 1
    set(f1,'PaperUnits','inches');
    set(f1, 'PaperPositionMode', 'manual');
    set(f1, 'PaperPosition', [-0.9 0 9 8]);
    set(f1,'PaperSize',[7.5 7.7]);   
    print(f1,'-dpdf','-r300','microring_mode_profiles.pdf');
end      