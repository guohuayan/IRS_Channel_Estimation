close all
clear all

load('elaa_covariance.mat','M','N','K','Cg','Cu_w','Cd_w','Loss_IU','Loss_BU');
% load('largescale_covariance.mat','M','N','K','Cg','Cu_w','Cd_w','Loss_IU','Loss_BU');
penetration=20;
Loss_BU=Loss_BU+penetration;
noise=-170+10*log10(200*1e3);
gain_BU=-Loss_BU-noise
gain_IU=-Loss_IU-noise
gain_BU=10.^(gain_BU./10);
gain_IU=10.^(gain_IU./10);
for k0=1:K
    Cd_w(:,:,k0)=Cd_w(:,:,k0).*gain_BU(k0);
    Cu_w(:,:,k0)=Cu_w(:,:,k0).*gain_IU(k0);
end
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
%         cond(C_IRS(:,:,m0,k0))
    end
end
iC_IRS0=zeros(N,N,M,K);
for m0=1:M
    for k0=1:K
        iC_IRS0(:,:,m0,k0)=inv(C_IRS(:,:,m0,k0));
    end
end
%% prepare pilot and phase configuration
tL=ceil((N-1)/K)+ceil(N/M);
% tL=N;
phi=dftmtx(N);
% theta0=ones(N,1);
% phi=diag(theta0)*phi;
pilot=dftmtx(K);
%%
% L1=ceil(N/M);
% gain=(L1*(K-1)+N)/N;
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
    iC_IRS=delta^2.*iC_IRS0;
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
        %% 
        phi_idx = randperm(N,tL);
        phi_t=phi(:,phi_idx);
        %% Phase II
        hat_Hd=Hd_w;
        Y2=zeros(M,K,N);
        Y2_SIC=zeros(M,K,tL);
        for l0=1:tL
            t_theta=phi_t(:,l0);
            tmp=(Hd_w+G.'*diag(t_theta)*hrw)*pilot+delta.*complex_randn(M,K);
            Y2(:,:,l0)=tmp/pilot;
            Y2_SIC(:,:,l0)=Y2(:,:,l0)-hat_Hd;
        end
        R=permute(Y2_SIC,[1,3,2]);
        %%
        LS_Hi=zeros(M,N,K);
        for k0=1:K
            R0=R(:,:,k0);
            bar_Rn=R0.';
            hat_Hg=zeros(N,M);
            Hmea=phi_t.';
            Cg_ref=C_IRS(:,:,:,k0);
            for m0=1:M
                Cgm=Cg_ref(:,:,m0);
                hat_Hg(:,m0)=Cgm*Hmea'/(Hmea*Cgm*Hmea'+delta^2/K*eye(tL))*bar_Rn(:,m0);
%                 hat_Hg(:,m0)=Cgm*Hmea'/(Hmea*Cgm*Hmea'+delta^2*eye(N)+var_Hd(m0,m0,1)*ones(N,N))*bar_Rn(:,m0);
            end
            LS_Hi(:,:,k0)=hat_Hg.';
        end
        Hg=zeros(M,N); %% M,K,tL
        for k0=1:K
            Hg=Hg+LS_Hi(:,:,k0);
        end
        Hg=Hg.';
        R=Y2_SIC;  %% M,K,tL
        Yr=zeros(tL*M,K);
        Gr=zeros(tL*M,N);
        for l0=1:tL
            t_theta=phi_t(:,l0);
            idx=(1:M)+(l0-1)*M;
            Yr(idx,:)=R(:,:,l0);
            Gr(idx,:)=Hg.'*diag(t_theta);
        end
        Hu=(Gr'*Gr)\Gr'*Yr;
        sigma_h=1/K;
        f0 = BALS_map_obj(R,Hg,Hu,phi_t,iC_IRS,sigma_h);
%         f0
        %%
        while(1)
            Hu_old=Hu;
            Hg_old=Hg;           
            Hu = BALS_Hu_map(R,Hg,phi_t,iC_IRS,sigma_h);
            Hg = BALS_Hg_map(R,Hu,phi_t,iC_IRS,sigma_h);
%             Hu = BALS_Hu(R,Hg,phi_t,Hu);
%             Hg = BALS_Hg(R,Hu,phi_t,Hg);
            f1 = BALS_map_obj(R,Hg,Hu,phi_t,iC_IRS,sigma_h);
%             f1
            if f1<f0
                Hu=Hu_old;
                Hg=Hg_old;
                fprintf('err\n');
                break
            elseif abs(log(abs(f1-f0)/abs(f0)))>5
                break
            end
            f0=f1;
        end
        MAP_Hi=zeros(M,N,K);
        for k0=1:K
            MAP_Hi(:,:,k0)=Hg.'*diag(Hu(:,k0));
        end
        %%
%         err=MAP_Hi-LS_Hi;
%         sum(abs(err(:)).^2)
        %%
        err=MAP_Hi-Hi_true;
        %%
        MSE_I_MAP_ite(i0)=sum(abs(err(:)).^2);
    end
%     NMSE_D(s0)=mean(MSE_D_ite)./mean(power_D_ite);
%     NMSE_I_LS(s0)=mean(MSE_I_LS_ite)./mean(power_I_ite);
    NMSE_I_MAP(s0)=mean(MSE_I_MAP_ite)./mean(power_I_ite);
end
% figure
% semilogy(snr_w,NMSE_D,'bo-');
% grid on
figure
% semilogy(snr_w,NMSE_I_LS,'bo-');
semilogy(snr_w,NMSE_I_MAP,'r+-');
grid on
save('base_BALS_map.mat','NMSE_I_MAP','snr_w');





