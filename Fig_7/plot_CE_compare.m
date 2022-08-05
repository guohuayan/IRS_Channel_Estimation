close all
clear all

load('baseline_improve.mat','NMSE_D','NMSE_I_LS','NMSE_I_MAP','snr_w');
% base_NMSE_D=NMSE_D;
% base_NMSE_I_LS=NMSE_I_LS;
basehd_NMSE_I_MAP=NMSE_I_MAP;
load('propose_cascade_final.mat','NMSE_D','NMSE_I_LS','NMSE_I_MAP','snr_w');
% pro_NMSE_D=NMSE_D;
% pro_NMSE_I_LS=NMSE_I_LS;
pro_NMSE_I_MAP=NMSE_I_MAP;
load('propose_cascade_optphi.mat','NMSE_D','NMSE_I_LS','NMSE_I_MAP','snr_w');
% opt_NMSE_I_LS=NMSE_I_LS;
opt_NMSE_I_MAP=NMSE_I_MAP;
%%
load('base_LMMSE_bound.mat','NMSE_I_MAP','snr_w');
base_LMMSE_MAP=NMSE_I_MAP;
load('base_BALS_bound.mat','NMSE_I_MAP','snr_w');
base_BALS_LS=NMSE_I_MAP;
load('base_BALS_map.mat','NMSE_I_MAP','snr_w');
base_BALS_MAP=NMSE_I_MAP;
load('baseline_improve_sameP.mat','NMSE_D','NMSE_I_LS','NMSE_I_MAP','snr_w');
baseimp_NMSE_I_MAP=NMSE_I_MAP;
%%
figure
semilogy(snr_w,base_BALS_LS,'m^-',snr_w,base_BALS_MAP,'ms-',snr_w,base_LMMSE_MAP,'mv-',...
    snr_w,basehd_NMSE_I_MAP,'kx-',snr_w,baseimp_NMSE_I_MAP,'k+-',...
    snr_w,pro_NMSE_I_MAP,'b+-',...
    snr_w,opt_NMSE_I_MAP,'ro-','LineWidth',1.1);
grid on
% ylim([1e-2,1e1]);
% xlim([10,30]);
xlabel('P_T (dBm)'); ylabel('NMSE');
% text(5,0.1,'Proposed (With phase-shift opt)','Color','r','FontSize',11);
% text(6,0.1,'Proposed (Without phase-shift opt)','Color','b','FontSize',11);
% text(20,0.8,'Baseline 5  (On-off with increased power)','Color','k','FontSize',11);
% text(5,0.8,'Baseline 4 (On-off)','Color','k','FontSize',11);
% text(20,0.8,'Baseline 1 (LMMSE)','Color','m','FontSize',11);
% text(20,2,'Baseline 3 (BALS-MAP)','Color','m','FontSize',11);
% text(20,7,'Baseline 2 (BALS)','Color','m','FontSize',11);
% x = [0.3,0.5];y = [0.6,0.5];
% annotation('arrow',x,y,'Color','r','Linewidth',1);
% x = [0.3,0.5];y = [0.3,0.4];
% annotation('arrow',x,y,'Color','b','Linewidth',1);
% x = [0.6,0.7];y = [0.3,0.4];
% annotation('arrow',x,y,'Color','k','Linewidth',1);
% x = [0.1,0.2];y = [0.3,0.4];
% annotation('arrow',x,y,'Color','k','Linewidth',1);
% x = [0.1,0.2];y = [0.5,0.6];
% annotation('arrow',x,y,'Color','m','Linewidth',1);
% x = [0.1,0.2];y = [0.5,0.6];
% annotation('arrow',x,y,'Color','m','Linewidth',1);
% x = [0.1,0.2];y = [0.7,0.8];
% annotation('arrow',x,y,'Color','m','Linewidth',1);
legend('Baseline 2 (BALS)','Baseline 3 (BALS-MAP)','Baseline 1 (LMMSE)',...
'Baseline 4 (On-off)','Baseline 5 (On-off with increased power)','Proposed (Without phase-shift opt)','Proposed (With phase-shift opt)',...
'location','southwest','FontSize',8);
