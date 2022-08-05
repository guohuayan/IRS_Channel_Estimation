close all
clear all

load('elaa_covariance.mat','M','N','K','Cg','Cu_w','Cd_w','Loss_IU','Loss_BU');
% load('largescale_covariance.mat','M','N','K','Cg','Cu_w','Cd_w','Loss_IU','Loss_BU');
penetration=20;
Loss_BU=Loss_BU+penetration;
noise=-170+10*log10(200*1e3);
gain_BU=-Loss_BU-noise;
gain_IU=-Loss_IU-noise;
gain_BU=10.^(gain_BU./10);
gain_IU=10.^(gain_IU./10);
for k0=1:K
    Cd_w(:,:,k0)=Cd_w(:,:,k0).*gain_BU(k0);
    Cu_w(:,:,k0)=Cu_w(:,:,k0).*gain_IU(k0);
end
%% fack Cg,Cu and calculate the inverse
TCg=zeros(N,N,M);
for m0=1:M
    row=(1:N)+(m0-1)*N;
    TCg(:,:,m0)=Cg(row,row);
end
TCu=Cu_w;
C_IRS=zeros(N,N,M,K); %% The covariance available
for m0=1:M
    for k0=1:K
        tmp=TCg(:,:,m0).*TCu(:,:,k0);
        tmp=tmp-diag(tmp)+real(diag(tmp));
        tmp=(tmp+tmp')/2;
        C_IRS(:,:,m0,k0)=tmp;
    end
end
iC_IRS0=zeros(N,N,M,K);
for m0=1:M
    for k0=1:K
        iC_IRS0(:,:,m0,k0)=inv(C_IRS(:,:,m0,k0));
    end
end
%% prepare pilot and phase configuration
phi_g=dftmtx(N);
pilot=dftmtx(K);
bar_x=pilot(:,1);
L1=ceil(N/M);
L2=N-L1;
% theta0=ones(N,1);
% phi=diag(theta0)*phi;
% phi1=phi(:,1:L1);
% phi2=phi(:,(L1+1):end);
%% simulation
ite=1e2;
snr_w=linspace(0,30,5);
NMSE_D=zeros(1,length(snr_w));
NMSE_D_LS=zeros(1,length(snr_w));
NMSE_I_LS=zeros(1,length(snr_w));
NMSE_I_MAP=zeros(1,length(snr_w));
% power_I=zeros(1,length(snr_w));
for s0=1:length(snr_w)
    rng(312);
    fprintf('s0=%d\n',s0);
    snr=snr_w(s0);
    delta=10^(-snr/20);
%     delta=1e-7;
    MSE_D_ite=zeros(1,ite);
    MSE_D_LS_ite=zeros(1,ite);
    MSE_I_LS_ite=zeros(1,ite);
    MSE_I_MAP_ite=zeros(1,ite);
    power_I_ite=zeros(1,ite);
    power_D_ite=zeros(1,ite);
    %%
%     iFCg=delta^2.*iFCg0;
%     iFCu=delta^2.*iFCu0;
    iC_IRS=delta^2.*iC_IRS0;
    srng = rng;
    rng(445);
    nrng=rng;
    for i0=1:ite  
        rng(srng);  
        %% generate channel
%         G=chol(Cg)*complex_randn(M*N,1);
        G=Cg^(1/2)*complex_randn(M*N,1);
        G=reshape(G,[N,M]);
        hrw=zeros(N,K);
        hd_w=zeros(M,K);
        for k0=1:K
%             hrw(:,k0)=chol(Cu_w(:,:,k0))*complex_randn(N,1);
            hrw(:,k0)=Cu_w(:,:,k0)^(1/2)*complex_randn(N,1);
%             hd_w(:,k0)=chol(Cd_w(:,:,k0))*complex_randn(M,1);
            hd_w(:,k0)=Cd_w(:,:,k0)^(1/2)*complex_randn(M,1);
        end
        Hi_true=zeros(M,N,K);
        for k0=1:K
            Hi_true(:,:,k0)=G.'*diag(hrw(:,k0));
        end
        srng = rng;
        %%
        rng(nrng);
%         rng('shuffle');
        theta0=exp(1j.*2*pi.*rand(N,1));
        phi=diag(theta0)*phi_g;
        phi1=phi(:,1:L1);
        phi2=phi(:,(L1+1):end);
%         rng(nrng);
        %% generate observation
        t_theta0=-phi1(:,1);
        Y0=(hd_w+G.'*diag(t_theta0)*hrw)*pilot+delta.*complex_randn(M,K);
        Y1_w=zeros(M,K,L1);
        for l0=1:L1
            t_theta=phi1(:,l0);
            Y1=(hd_w+G.'*diag(t_theta)*hrw)*pilot+delta.*complex_randn(M,K);
            Y1_w(:,:,l0)=Y1;
        end
        Y2=zeros(M,L2);
        for l0=1:L2
            t_theta=phi2(:,l0);
            y2=(hd_w+G.'*diag(t_theta)*hrw)*bar_x+delta.*complex_randn(M,1);
            Y2(:,l0)=y2;
        end
        nrng = rng;
        %% estimate direct channel
        R0=(Y0+Y1_w(:,:,1))/2/pilot;
        hat_Hd_w=zeros(M,K);
        for k0=1:K
            Cdk=Cd_w(:,:,k0);
            hat_Hd_w(:,k0)=Cdk/(Cdk+delta^2/2/K*eye(M))*R0(:,k0);
        end
        err=hat_Hd_w-hd_w;
        power_d=sum(abs(hd_w(:)).^2);
        power_D_ite(i0)=power_d;
%         MSE_D_ite(i0)=sum(abs(err(:)).^2)./power_d;
        MSE_D_ite(i0)=sum(abs(err(:)).^2);
%         sum(abs(err(:)).^2)
        err=R0-hd_w;
%         MSE_D_LS_ite(i0)=sum(abs(err(:)).^2)./power_d;
        MSE_D_LS_ite(i0)=sum(abs(err(:)).^2);
%         sum(abs(err(:)).^2)
        %% pre-processing IRS
        RL=zeros(M,K,L1);
        tRL=zeros(M,K,L1);
        for idx=1:L1
            if idx==1
                RL(:,:,idx)=(Y1_w(:,:,1)-Y0)/2;
            else
                RL(:,:,idx)=Y1_w(:,:,idx)-(Y0+Y1_w(:,:,1))/2;
            end
            tRL(:,:,idx)=RL(:,:,idx)/pilot;
        end
        sigma1_sq=1/2/K;
        sigmaL_sq=1/2/K*3;
        bar_r=zeros(M,L2);
        tmp=(Y0+Y1_w(:,:,1))/2;
        tmp=tmp(:,1);
        for l0=1:L2
            bar_r(:,l0)=Y2(:,l0)-tmp;
        end
        bar_sigmaL_sq=3/2;
        %% initial
        bar_R=zeros(M,N);
        for l0=1:L1
            tmp=RL(:,:,l0);
            bar_R(:,l0)=tmp(:,1);
        end
        bar_R(:,(L1+1):end)=bar_r;
        LS_Hg=bar_R/phi;
        measure=zeros(M*L1,N);
        observation=zeros(M*L1,K);
        for l0=1:L1
            idx0=(1:M)+M*(l0-1);
            t_theta=phi1(:,l0);
            measure(idx0,:)=LS_Hg*diag(t_theta);
            observation(idx0,:)=RL(:,:,l0);
        end
        LS_Hu=(measure'*measure)\measure'*observation/pilot;
        LS_Hi=zeros(M,N,K);
        for k0=1:K
            LS_Hi(:,:,k0)=LS_Hg*diag(LS_Hu(:,k0));
        end
        err=LS_Hi-Hi_true;
        power_i=sum(abs(Hi_true(:)).^2);
        MSE_I_LS_ite(i0)=sum(abs(err(:)).^2);
        power_I_ite(i0)=power_i;
        %%
        Hg=LS_Hg.';
        Hu=LS_Hu;
        [f0,~,~] = obj_fun_MAP_cascade(Hu,Hg,tRL,bar_r,bar_x,phi,iC_IRS,sigma1_sq,sigmaL_sq,bar_sigmaL_sq);
%         f0
        %%
        while(1)
            Hu_old=Hu;Hg_old=Hg;
            Hu = MAP_hu_optimal(Hg,tRL,bar_r,bar_x,phi,iC_IRS,sigma1_sq,sigmaL_sq,bar_sigmaL_sq);
            Hg = MAP_hg_propose_cascade(Hu,Hg,tRL,bar_r,bar_x,phi,iC_IRS,sigma1_sq,sigmaL_sq,bar_sigmaL_sq);
            [f1,~,~] = obj_fun_MAP_cascade(Hu,Hg,tRL,bar_r,bar_x,phi,iC_IRS,sigma1_sq,sigmaL_sq,bar_sigmaL_sq);
%             f1
            if f1<f0
                fprintf('err\n');
                Hu=Hu_old;
                Hg=Hg_old;
                break
            elseif abs(log(abs(f1-f0)/abs(f0)))>3
                break
            end
            f0=f1;
        end
        %%
%         pause
        MAP_Hi=zeros(M,N,K);
        for k0=1:K
            MAP_Hi(:,:,k0)=Hg.'*diag(Hu(:,k0));
        end
        err=MAP_Hi-Hi_true;
        MSE_I_MAP_ite(i0)=sum(abs(err(:)).^2);
    end
    NMSE_D(s0)=mean(MSE_D_ite)./mean(power_D_ite);
    NMSE_D_LS(s0)=mean(MSE_D_LS_ite)./mean(power_D_ite);
    NMSE_I_LS(s0)=mean(MSE_I_LS_ite)./mean(power_I_ite);
    NMSE_I_MAP(s0)=mean(MSE_I_MAP_ite)./mean(power_I_ite);
%     NMSE_I_LS(s0)=mean(MSE_I_LS_ite);
%     power_I=mean(power_I_ite);
end
% figure
% semilogy(snr_w,NMSE_D,'ro-',snr_w,NMSE_D_LS,'bo-');
% grid on
figure
semilogy(snr_w,NMSE_I_LS,'bo-',snr_w,NMSE_I_MAP,'r+-');
grid on
% figure
% semilogy(snr_w,power_I,'ro-');
% grid on
save('propose_cascade_final.mat','NMSE_D','NMSE_I_LS','NMSE_I_MAP','snr_w');
%%



