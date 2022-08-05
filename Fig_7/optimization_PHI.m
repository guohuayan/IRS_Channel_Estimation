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
phi=dftmtx(N);
L1=ceil(N/M);
phi1=phi(:,1:L1);
% theta0=ones(N,1);
% phi=diag(theta0)*phi;
%%
E=zeros(N,N);
for l0=1:L1
    t_phi=phi1(:,l0);
    for m0=1:M
        for k0=1:K
            Ctmp=diag(t_phi)'*conj(C_IRS(:,:,m0,k0))*diag(t_phi);
            E=E+Ctmp;
        end
    end
end
theta=ones(N,1);
f0=real(theta'*E*theta)
length_thr=20;
while(1)
    theta_old=theta;
    tmp=E*theta;
    theta=exp(1j.*angle(tmp));
    f1=real(theta'*E*theta);
    if f1<f0
        theta=theta_old;
        fprintf('err HU\n');
        break
    elseif abs(log(abs(f1-f0)/abs(f0)))>length_thr  %abs(f1-f0)<1e-6
        break
    end
    f0=f1;
end
f1
save('largescale_covariance_opt_phi.mat','M','N','K','Cg','Cu_w','Cd_w','C_IRS','iC_IRS0','theta');
%%




