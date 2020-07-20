% This  script uses a the 'axes' command instead of 'subplot' to create
% better figures.


close all;
clear all;
clc;

fs_xylabel = 25;
fs_ticklabel = 20;
fs_legend = 20;

x0 = 0.135;
y0 = 0.11;
right_clearence = 0.09;
top_clearence = 0.05;
hspacing = 0.075;
vspacing = 0.09;
W = (1 - x0 - hspacing - right_clearence)/2;
H = (1 - y0 - vspacing - top_clearence)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1: Resonant wavelength and Q versus m-number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


m_arr = 25:1:33;

lambda_TE_anal = [1769.29 1710.52 1655.62 1604.21 1555.95 1510.57 1467.81 1427.86 1389.3];

lambda_TE_FD = [1775.27 1716.25 1661.1 1609.47 1561.00 1515.436 1472.50 1431.97 1393.66 ];

lambda_TM_anal = [ 1825.33 1762.79 1704.49 1650 1598.96 1551.04 1505.96 1463.48 1423.37];

lambda_TM_FD = [1823.64 1761.13 1702.85 1648.38 1597.359 1549.461 1504.41 1461.947 1421.86];

Q_TE_anal = [4.31e11 1.34e12 4.21e12 1.32e13 4.14e13 1.30e14 4.11e14 1.29e15 4.11e15];

Q_TE_FD = [4.499e11 1.385e12 4.29e12 1.32e13 4.099e13 1.2766e14 3.963e14 1.23e15 3.84e15];

Q_TM_anal = [4.33e11 1.348e12 4.20e12 1.314e13 4.12e13 1.29248e14 4.07e14 1.21196e15 4.029e15];

Q_TM_FD = [4.643e11 1.462e12 4.549e12 1.432e13 4.5195e13 1.430836e14 4.5521e14 1.431e15 4.534e15];

%--------------------------------------------------------------------------

f1 = figure(1) ; 

f11 = axes( 'pos', [x0  y0 + H + vspacing W H]);

plot(m_arr,lambda_TE_anal,'b-',...
     m_arr,lambda_TE_FD,'ro',...
     'markerfacecolor','r',...
     'markeredgecolor','r',...
     'markersize',12,...
     'linewidth',1.5);

xlabel('$m$ number','Fontsize',fs_xylabel,'interpreter','latex');

ylabel('Wavelength (nm)','Fontsize',fs_xylabel,'interpreter','latex');

xlim([23 35]);

set(f11,'fontsize',fs_ticklabel);

text(23.2, 1340,'(a) TE Mode','fontsize',25,'fontname','helvetica');

legend('$\mathsf{Analytical}$',...
       '$\mathsf{Finite \, Difference}$');
set(legend(f11),'FontSize',fs_legend,...
                'EdgeColor',[0.99 0.99 0.99],...
                'Position',[0.30 0.83 0.2 0.13],...
                'interpreter','latex',...
                'Color','w');


%--------------------------------------------------------------------------

f12 = axes( 'pos',[x0 + W + hspacing y0 + H + vspacing W H]);

plot(m_arr,lambda_TM_anal,'b-',...
     m_arr,lambda_TM_FD,'ro',...
     'markerfacecolor','r',...
     'markeredgecolor','r',...
     'markersize',12,...
     'linewidth',1.5);
 
xlabel('$m$ number','Fontsize',fs_xylabel,'interpreter','latex');

ylabel('Wavelength (nm)','Fontsize',fs_xylabel,'interpreter','latex');

xlim([23 35]);

set(f12,'fontsize',fs_ticklabel);

text(23.2, 1440,'(c) TM Mode','fontsize',25,'fontname','helvetica');

legend('$\mathsf{Analytical}$',...
       '$\mathsf{Finite \, Difference}$');
set(legend(f12),'FontSize',fs_legend,...
                'EdgeColor',[0.99 0.99 0.99],...
                'Position',[0.725 0.83 0.2 0.13],...
                'interpreter','latex',...
                'Color','w');
   

%--------------------------------------------------------------------------

f13 = axes('pos', [x0  y0 W H]);

semilogy(m_arr,Q_TE_anal,'b-',...
         m_arr,Q_TE_FD,'ro',...
         'markerfacecolor','r',...
         'markeredgecolor','r',...
         'markersize',12,...
         'linewidth',1.5); 

xlabel('$m$ number','Fontsize',fs_xylabel,'interpreter','latex');

ylabel('$Q$','Fontsize',fs_xylabel,'interpreter','latex');

xlim([23 35]);

set(f13,'fontsize',fs_ticklabel);

text(23.2, 4e15,'(b) TE Mode','fontsize',25,'fontname','helvetica');

legend('$\mathsf{Analytical}$',...
       '$\mathsf{Finite \, Difference}$');
set(legend(f13),'FontSize',fs_legend,...
                'EdgeColor',[0.99 0.99 0.99],...
                'Position',[0.30 0.095 0.2 0.13],...
                'interpreter','latex',...
                'Color','w');

%--------------------------------------------------------------------------

f14 = axes( 'pos',[x0 + W + hspacing  y0 W H]);

semilogy(m_arr,Q_TM_anal,'b-',...
         m_arr,Q_TM_FD,'ro',...
         'markerfacecolor','r',...
         'markeredgecolor','r',...
         'markersize',12,...
         'linewidth',1.5);

xlabel('$m$ number','Fontsize',fs_xylabel,'interpreter','latex');

ylabel('$Q$','Fontsize',fs_xylabel,'interpreter','latex');

xlim([23 35]);

set(f14,'fontsize',fs_ticklabel);

text(23.2, 4e15,'(d) TE Mode','fontsize',25,'fontname','helvetica');

legend('$\mathsf{Analytical}$',...
       '$\mathsf{Finite \, Difference}$');
set(legend(f14),'FontSize',fs_legend,...
                'EdgeColor',[0.99 0.99 0.99],...
                'Position',[0.725 0.095 0.2 0.13],...
                'interpreter','latex',...
                'Color','w');


print_pdf_1 = 1;

if print_pdf_1 == 1
    set(f1,'PaperUnits','inches');
    set(f1, 'PaperPositionMode', 'manual');
    set(f1, 'PaperPosition', [-1.2 -0.2 17.2 9.0]);
    set(f1,'PaperSize',[14.5 8.5]);   
    print(f1,'-dpdf','-r300','FD_vs_analytic_mvar.pdf');
    %print(f1,'-deps','-r300','FD_vs_analytic_mvar.eps');
end


%--------------------------------------------------------------------------
%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2: Resonant wavelength and Q versus grid-size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gridsize_arr = [200 150 100 50 25 20];

lambda_TE_FD_gs = [1635.901 1603.037 1570.5658 1565.15 1562.06 1561.00];
lambda_TE_anal_gs = 1555.95 * ones(size(gridsize_arr));

Q_TE_gs = [3.55e12 6.50e12 1.23e13 2.66e13 3.83e13 4.099e13];
 
lambda_TM_FD_gs = [1643.153 1614.1219 1588.293 1593.22 1596.931 1597.359];

Q_TM_gs = [2.73e13 3.395e13 4.2e13 4.72e13 4.60e13 4.52e13];

f2 = figure(2); 
%--------------------------------------------------------------------------

f21 = axes( 'pos', [x0  y0 + H + vspacing W H]);

semilogx(gridsize_arr,lambda_TE_FD_gs,'o',...
      'markerfacecolor',[0.5 0 0],...
      'markeredgecolor',[0.5 0 0],...
      'markersize',12); hold on;

line([10 300],[1555.95 1555.95],'LineStyle','--','linewidth',1.5);

xlim([12 250]);
ylim([1539 1650]);

set(f21,'xdir','reverse',...
     'ydir','reverse',...
     'XTick',[10 20 25 50 100 200],...
     'XTicklabel','10|20|25|50|100|200',...
     'fontsize',fs_ticklabel);

xlabel('Grid size (nm)','interpreter','latex','fontsize',fs_xylabel);
ylabel('Wavelength (nm)','interpreter','latex','fontsize',fs_xylabel);

text(34, 1640,'(a) TE Mode','fontsize',25,'fontname','helvetica');
text(240, 1550,'Analytical','fontsize',20,'fontname','helvetica');

%--------------------------------------------------------------------------

f22 = axes( 'pos',[x0 + W + hspacing y0 + H + vspacing W H]);

semilogx(gridsize_arr,lambda_TM_FD_gs,'o',...
      'markerfacecolor',[0.5 0 0],...
      'markeredgecolor',[0.5 0 0],...
      'markersize',12); hold on;

line([10 300],[1598.96 1598.96],'LineStyle','--','linewidth',1.5);

xlim([12 250]);
ylim([1539 1650]);

set(gca,'xdir','reverse',...
     'ydir','reverse',...
     'XTick',[10 20 25 50 100 200],...
     'XTicklabel','10|20|25|50|100|200',...
     'fontsize',fs_ticklabel);

xlabel('Grid size (nm)','interpreter','latex','fontsize',fs_xylabel);
ylabel('Wavelength (nm)','interpreter','latex','fontsize',fs_xylabel);

text(34, 1640,'(c) TM Mode','fontsize',25,'fontname','helvetica');
text(240, 1593,'Analytical','fontsize',20,'fontname','helvetica');

%--------------------------------------------------------------------------

f23 = axes('pos', [x0  y0 W H]);

loglog(gridsize_arr,Q_TE_gs,'o',...
      'markerfacecolor',[0.5 0 0],...
      'markeredgecolor',[0.5 0 0],...
      'markersize',12); hold on;

line([10 300],[4.14e13 4.14e13],'LineStyle','--','linewidth',1.5);

xlim([12 250]);
ylim([1e12 1e14]);

set(gca,'xdir','reverse',...
     'XTick',[10 20 25 50 100 200],...
     'XTicklabel','10|20|25|50|100|200',...
     'fontsize',fs_ticklabel);

xlabel('Grid size (nm)','interpreter','latex','fontsize',fs_xylabel);
ylabel('$Q$','interpreter','latex','fontsize',fs_xylabel);

text(34, 1.5e12,'(b) TE Mode','fontsize',25,'fontname','helvetica');
text(240, 5.3e13,'Analytical','fontsize',20,'fontname','helvetica');

%--------------------------------------------------------------------------

f24 = axes( 'pos',[x0 + W + hspacing  y0 W H]);
loglog(gridsize_arr,Q_TM_gs,'o',...
      'markerfacecolor',[0.5 0 0],...
      'markeredgecolor',[0.5 0 0],...
      'markersize',12); hold on;

line([10 300],[4.12e13 4.12e13],'LineStyle','--','linewidth',1.5);

xlim([12 250]);
ylim([1e12 1e14]);

set(gca,'xdir','reverse',...
        'XTick',[10 20 25 50 100 200],...
        'XTicklabel','10|20|25|50|100|200',...
        'fontsize',fs_ticklabel);

xlabel('Grid size (nm)','interpreter','latex','fontsize',fs_xylabel);
ylabel('$Q$','interpreter','latex','fontsize',fs_xylabel);

text(34, 1.5e12,'(d) TM Mode','fontsize',25,'fontname','helvetica');
text(240, 5.3e13,'Analytical','fontsize',20,'fontname','helvetica');

print_pdf_2 = 1;

if print_pdf_2 == 1
    set(f2,'PaperUnits','inches');
    set(f2, 'PaperPositionMode', 'manual');
    set(f2, 'PaperPosition', [-1.2 -0.2 17.2 9.0]);
    set(f2,'PaperSize',[14.5 8.5]);   
    print(f2,'-dpdf','-r300','FD_vs_analytic_gridvar.pdf');
    %print(f2,'-deps','-r300','FD_vs_analytic_gridvar.eps');
end

