close all
clear all

load('elaa_covariance.mat','M','N','K','Cg','Cu_w','Cd_w','Loss_IU','Loss_BU');
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
%% simulation
ite=1e2;
snr=15;
loop=10;
NMSE_I_MAP=zeros(loop,ite);
objfun=zeros(loop,ite);
rng(312);
delta=10^(-snr/20);
power_I_ite=zeros(1,ite);
iC_IRS=delta^2.*iC_IRS0;
srng = rng;
rng(445);
nrng=rng;
for i0=1:ite  
    if mod(i0,20)==0
        fprintf('ite=%d\n',i0);
    end
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
    power_i=sum(abs(Hi_true(:)).^2);
    power_I_ite(i0)=power_i;
    srng = rng;
    rng(nrng);
    %%
    theta0=exp(1j.*2*pi.*rand(N,1));
    phi=diag(theta0)*phi_g;
    phi1=phi(:,1:L1);
    phi2=phi(:,(L1+1):end);
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
    %%
    Hg=LS_Hg.';
    Hu=LS_Hu;
    %%
    [f0,~,~] = obj_fun_MAP_cascade(Hu,Hg,tRL,bar_r,bar_x,phi,iC_IRS,sigma1_sq,sigmaL_sq,bar_sigmaL_sq);
%     MAP_Hi=zeros(M,N,K);
%     for k0=1:K
%         MAP_Hi(:,:,k0)=Hg.'*diag(Hu(:,k0));
%     end
%     err=MAP_Hi-Hi_true;
%     NMSE_I_MAP(1,i0)=sum(abs(err(:)).^2);
%     objfun(1,i0)=f0;
        %%
    for lp=1:loop
%             Hu_old=Hu;Hg_old=Hg;
        Hu = MAP_hu_optimal(Hg,tRL,bar_r,bar_x,phi,iC_IRS,sigma1_sq,sigmaL_sq,bar_sigmaL_sq);
        [f0,~,~] = obj_fun_MAP_cascade(Hu,Hg,tRL,bar_r,bar_x,phi,iC_IRS,sigma1_sq,sigmaL_sq,bar_sigmaL_sq);
        objfun(lp,i0)=f0;
        Hg = MAP_hg_propose_cascade(Hu,Hg,tRL,bar_r,bar_x,phi,iC_IRS,sigma1_sq,sigmaL_sq,bar_sigmaL_sq);
%         MAP_Hi=zeros(M,N,K);
%         for k0=1:K
%             MAP_Hi(:,:,k0)=Hg.'*diag(Hu(:,k0));
%         end
%         err=MAP_Hi-Hi_true;
%         NMSE_I_MAP(lp,i0)=sum(abs(err(:)).^2);
%         Hg = MAP_hg_propose_cascade(Hu,Hg,tRL,bar_r,bar_x,phi,iC_IRS,sigma1_sq,sigmaL_sq,bar_sigmaL_sq);
    end
end
% NMSE_I_MAP=mean(NMSE_I_MAP,2);
% power=mean(power_I_ite);
% NMSE_I_MAP=NMSE_I_MAP/power;
objfun=mean(objfun,2);
figure
% semilogy(1:loop,NMSE_I_MAP,'r+-');
semilogy(1:loop,objfun,'r+-');
grid on
% save('iteration_final_30.mat','objfun','loop');
load('iteration_final_15_n.mat','objfun','loop');
%%



