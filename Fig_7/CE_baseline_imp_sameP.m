%% increase the power such that the power is equal to the proposed algorithm
close all
clear all

% load('elaa_covariance.mat','M','N','K','Cg','Cu_w','Cd_w','Loss_IU','Loss_BU');
% penetration=20;
% Loss_BU=Loss_BU+penetration;
% noise=-170+10*log10(200*1e3);
% gain_BU=-Loss_BU-noise;
% gain_IU=-Loss_IU-noise;
% gain_BU=10.^(gain_BU./10);
% gain_IU=10.^(gain_IU./10);
% for k0=1:K
%     Cd_w(:,:,k0)=Cd_w(:,:,k0).*gain_BU(k0);
%     Cu_w(:,:,k0)=Cu_w(:,:,k0).*gain_IU(k0);
% end
load('baseline_elaa_covariance.mat','M','N','K','Cg','Cu_w','Cd_w','C_lambda');
TCg=zeros(N,N,M);
for m0=1:M
    row=(1:N)+(m0-1)*N;
    TCg(:,:,m0)=Cg(row,row);
end
TCu=Cu_w;
C_IRS=zeros(N,N,M,K); %% The covariance available
for m0=1:M
    for k0=1:K
        C_IRS(:,:,m0,k0)=TCg(:,:,m0).*TCu(:,:,k0);
    end
end
%% prepare pilot and phase configuration
phi=dftmtx(N);
theta0=ones(N,1);
phi=diag(theta0)*phi;
pilot=dftmtx(K);
%%
plen_1=N;
plen_2=ceil(N/M)*(K-1);
plen=plen_1+plen_2;
p_new_1=plen/plen_1;
p_new_2=plen/plen_2;
%% simulation
ite=1e2;
snr_w=linspace(0,30,5);
NMSE_D=zeros(1,length(snr_w));
NMSE_D_LS=zeros(1,length(snr_w));
NMSE_I_LS=zeros(1,length(snr_w));
NMSE_I_MAP=zeros(1,length(snr_w));
for s0=1:length(snr_w)
    rng(312);
    fprintf('s0=%d\n',s0);
    snr=snr_w(s0);
    delta=10^(-snr/20);
%     delta=1e-7;
    MSE_D_ite=zeros(1,ite);
    MSE_I_LS_ite=zeros(1,ite);
    MSE_I_MAP_ite=zeros(1,ite);
    power_I_ite=zeros(1,ite);
    power_D_ite=zeros(1,ite);
    %%
    srng = rng;
    for i0=1:ite  
        rng(srng); 
        %% generate channel
%         G=chol(Cg)*complex_randn(M*N,1);
        G=Cg^(1/2)*complex_randn(M*N,1);
        G=reshape(G,[N,M]);
        hrw=zeros(N,K);
        Hd_w=zeros(M,K);
        for k0=1:K
            hrw(:,k0)=Cu_w(:,:,k0)^(1/2)*complex_randn(N,1);
            Hd_w(:,k0)=Cd_w(:,:,k0)^(1/2)*complex_randn(M,1);
        end
        Hi_true=zeros(M,N,K);
        for k0=1:K
            Hi_true(:,:,k0)=G.'*diag(hrw(:,k0));
        end
        srng = rng;
        power_d=sum(abs(Hd_w(:)).^2);
        power_D_ite(i0)=power_d;
        power_i=sum(abs(Hi_true(:)).^2);
        power_I_ite(i0)=power_i;
        %% Phase I
        Y1=Hd_w*pilot+delta.*complex_randn(M,K);
        R0=Y1/pilot;
        hat_Hd=zeros(M,K);
        var_Hd=zeros(M,M,K);
        for k0=1:K
            Cdk=Cd_w(:,:,k0);
            hat_Hd(:,k0)=Cdk/(Cdk+delta^2/K*eye(M))*R0(:,k0);
            var_Hd(:,:,k0)=inv(inv(Cdk)+eye(M)/(delta^2/K));
        end
        err=hat_Hd-Hd_w;
        MSE_D_ite(i0)=sum(abs(err(:)).^2);
        %% Phase II
        G=G.';
        a=ones(1,N);
        Y2=Hd_w(:,1)*a+G*diag(hrw(:,1))*phi+delta.*complex_randn(M,N)./sqrt(p_new_1);
        Y2_SIC=Y2-hat_Hd(:,1)*a;
        LS_G=Y2_SIC/phi;
        %%
        Y2_SIC=Y2-hat_Hd(:,1)*a;
        bar_Rn=Y2_SIC.';
        hat_Hg=zeros(N,M);
        Hmea=phi.';
        Cg_ref=C_IRS(:,:,:,1);
        for m0=1:M
            Cgm=Cg_ref(:,:,m0);
            hat_Hg(:,m0)=Cgm*Hmea'/(Hmea*Cgm*Hmea'+delta^2*eye(N)/p_new_1+var_Hd(m0,m0,1)*(a*a'))*bar_Rn(:,m0);
%             hat_Hg(:,m0)=Cgm*Hmea'/(Hmea*Cgm*Hmea'+delta^2*eye(N))*bar_Rn(:,m0);
        end
        LS_G=hat_Hg.';
        %% phase_III
%         hat_Hu=zeros(N,K);
%         hat_Hu(:,1)=a.';
%         tau=N/M;
%         for k0=2:K
%             hatHu_k=zeros(N,1);
%             for n0=1:tau
%                 active=(1:M).'+(n0-1)*M;
%                 Y_tmp=Hd_w(:,k0)+G(:,active)*hrw(active,k0)+delta.*complex_randn(M,1)./sqrt(K-1)./sqrt(p_new_2);
%                 Y_tmp=Y_tmp-hat_Hd(:,k0);
%                 hatHu_k(active)=(LS_G(:,active)'*LS_G(:,active))\LS_G(:,active)'*Y_tmp;
%             end
%             hat_Hu(:,k0)=hatHu_k;
%         end
%         LS_Hi=zeros(M,N,K);
%         LS_Hi(:,:,1)=LS_G;
%         for k0=2:K
%             LS_Hi(:,:,k0)=LS_G*diag(hat_Hu(:,k0));
%         end
%         err=LS_Hi-Hi_true;
%         MSE_I_LS_ite(i0)=sum(abs(err(:)).^2);
        %%
% %         hat_Hd=Hd_w;
%         Y2_SIC=Y2-hat_Hd(:,1)*a;
%         bar_Rn=Y2_SIC.';
%         hat_Hg=zeros(N,M);
%         Hmea=phi.';
%         Cg_ref=C_IRS(:,:,:,1);
%         for m0=1:M
%             Cgm=Cg_ref(:,:,m0);
%             hat_Hg(:,m0)=Cgm*Hmea'/(Hmea*Cgm*Hmea'+delta^2*eye(N)+var_Hd(m0,m0,1)*(a*a'))*bar_Rn(:,m0);
% %             hat_Hg(:,m0)=Cgm*Hmea'/(Hmea*Cgm*Hmea'+delta^2*eye(N))*bar_Rn(:,m0);
%         end
%         LS_G=hat_Hg.';
        hat_Hu=zeros(N,K);
        hat_Hu(:,1)=a.';
        tau=N/M;
        for k0=2:K
            hatHu_k=zeros(N,1);
            C_hu=C_lambda(:,:,k0-1);
            for n0=1:tau
                active=(1:M).'+(n0-1)*M;
                Y_tmp=Hd_w(:,k0)+G(:,active)*hrw(active,k0)+delta.*complex_randn(M,1)./sqrt(K-1)./sqrt(p_new_2);
                Y_tmp=Y_tmp-hat_Hd(:,k0);
                Hmea=LS_G(:,active);
                Cgm=C_hu(active,active);
                hatHu_k(active)=Cgm*Hmea'/(Hmea*Cgm*Hmea'+delta^2*eye(M)/(K-1)/p_new_2+var_Hd(:,:,k0))*Y_tmp;
%                 hatHu_k(active)=Cgm*Hmea'/(Hmea*Cgm*Hmea'+delta^2*eye(M))*Y_tmp;
            end
            hat_Hu(:,k0)=hatHu_k;
        end
        LS_Hi=zeros(M,N,K);
        LS_Hi(:,:,1)=LS_G;
        for k0=2:K
            LS_Hi(:,:,k0)=LS_G*diag(hat_Hu(:,k0));
        end
        err=LS_Hi-Hi_true;
        MSE_I_MAP_ite(i0)=sum(abs(err(:)).^2);
    end
    NMSE_D(s0)=mean(MSE_D_ite)./mean(power_D_ite);
    NMSE_I_LS(s0)=mean(MSE_I_LS_ite)./mean(power_I_ite);
    NMSE_I_MAP(s0)=mean(MSE_I_MAP_ite)./mean(power_I_ite);
end
figure
semilogy(snr_w,NMSE_D,'bo-');
grid on
figure
semilogy(snr_w,NMSE_I_LS,'bo-',snr_w,NMSE_I_MAP,'r+-');
grid on
save('baseline_improve_sameP.mat','NMSE_D','NMSE_I_LS','NMSE_I_MAP','snr_w');





