close all
clear all

load('iteration_final_15_n.mat','objfun','loop');
delta=10^(15/10);
final_obj_15=objfun*delta;
load('iteration_opt_15_n.mat','objfun','loop');
opt_obj_15=objfun*delta;
figure
semilogy(1:loop,final_obj_15,'b+-',1:loop,opt_obj_15,'ro-',...
    'LineWidth',1.1);
grid on
xlabel('Iterations'); ylabel('Objective Function');
% text(5,-1e5,'Proposed (With phase-shift opt)','Color','r','FontSize',11);
% text(6,-1e4,'Proposed (Without phase-shift opt)','Color','b','FontSize',11);
% x = [0.3,0.5];y = [0.6,0.5];
% annotation('arrow',x,y,'Color','r','Linewidth',1);
% x = [0.3,0.5];y = [0.3,0.4];
% annotation('arrow',x,y,'Color','b','Linewidth',1);
legend('Proposed (Without phase-shift opt)','Proposed (With phase-shift opt)','FontSize',11,'location','southeast');