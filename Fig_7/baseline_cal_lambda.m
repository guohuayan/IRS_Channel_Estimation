close all
clear all

load('elaa_covariance.mat','M','N','K','Cg','Cu_w','Cd_w','Loss_IU','Loss_BU');
% load('largescale_covariance.mat','M','N','K','Cg','Cu_w','Cd_w','Loss_IU','Loss_BU');
% penetration=20;
% Loss_BU=Loss_BU+penetration;
noise=-170+10*log10(200*1e3);
% gain_BU=-Loss_BU-noise;
gain_IU=-Loss_IU-noise;
% gain_BU=10.^(gain_BU./10);
gain_IU=10.^(gain_IU./10);
for k0=1:K
%     Cd_w(:,:,k0)=Cd_w(:,:,k0).*gain_BU(k0);
    Cu_w(:,:,k0)=Cu_w(:,:,k0).*gain_IU(k0);
end
%% prepare pilot and phase configuration
ite=1e4;
C_lambda=zeros(N,N,K-1);
for i0=1:ite  
    %% generate channel
    hrw=zeros(N,K);
    lambda=zeros(N,K-1);
    hrw(:,1)=Cu_w(:,:,1)^(1/2)*complex_randn(N,1);
    for k0=2:K
        idx=k0-1;
        hrw(:,k0)=Cu_w(:,:,k0)^(1/2)*complex_randn(N,1);
        lambda(:,idx)=hrw(:,k0)./hrw(:,1);
        C_lambda(:,:,idx)=C_lambda(:,:,idx)+lambda(:,idx)*lambda(:,idx)';
    end
end
C_lambda=C_lambda./ite;
save('baseline_elaa_covariance.mat','N','K','C_lambda');