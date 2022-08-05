close all
clear all

rng(123);
%% 
% fc=5.25;
fc=2.4;
lightspeed=3e8;
%% point to point
DS_mu=(-0.28*log10(1+fc)-7.173);
DS_std=(0.1*log10(1+fc)+0.055);
ASD_mu=1.62;
ASD_std=0.25;
ASA_mu=(-0.11*log10(1+fc)+1.863);
ASA_std=(0.12*log10(1+fc)+0.059);
ZSA_mu=(-0.15*log10(1+fc)+1.387);
ZSA_std=(-0.09*log10(1+fc)+0.746);
ZSD_mu=1.08;
ZSD_std=0.36;
%%
c_ASD=5;
c_ASA=11;
c_ZSA=9;
c_ZSD=3/8*10^(ZSD_mu);
%%
SFc=3;
Nray=30;
Thr=10^(-25/10);
% Thr=10^(-30/10);
%% IRS Arrive
IAS_mu=ASA_mu;
IAS_std=ASA_std;
IZS_mu=ZSA_mu;
IZS_std=ZSA_std;
c_A_IRS=c_ASA;
c_Z_IRS=c_ZSA;
%%
% Ncd=40; 
Ncd=19; 
DS=10^DS_mu;
ASD=min(10^ASD_mu,104);
ZSD=min(10^ZSD_mu,52);
ASI=min(10^IAS_mu,104);
ZSI=min(10^IZS_mu,52);
%% BS-to-IRS
% while(1)
    tau_p=rand(Ncd,1).*2.*10^(DS_mu+DS_std);
    tau=tau_p-min(tau_p);
    phi_AoD=2*10^(ASD_mu+ASD_std).*(1-2.*rand(Ncd,1));
    theta_ZoD=2*10^(ZSD_mu+ZSD_std).*(1-2.*rand(Ncd,1));
    phi_AoA=2*10^(IAS_mu+IAS_std).*(1-2.*rand(Ncd,1));
    theta_ZoA=2*10^(IZS_mu+IZS_std).*(1-2.*rand(Ncd,1));
    Zn=SFc.*randn(Ncd,1);
    Pn_p=exp(-tau./DS).*10.^(-Zn./10).*exp(-abs(phi_AoD./ASD)*sqrt(2)).*exp(-abs(theta_ZoD./ZSD)*sqrt(2))...
        .*exp(-abs(phi_AoA./ASI)*sqrt(2)).*exp(-abs(theta_ZoA./ZSI)*sqrt(2));
%     Pn_p=exp(-tau./DS).*10.^(-Zn./10).*exp(-abs(phi_AoD./ASD).*sqrt(2))...
%         .*exp(-abs(phi_AoA./ASI)*sqrt(2)).*exp(-abs(theta_ZoA./ZSI)*sqrt(2));
    Pn=Pn_p./sum(Pn_p);
%     [Pn, Ncdt] = cluster_elimination(Pn,Thr);
%     Ncdt
%     if Ncdt>20
%         break
%     end
% end
% figure
% semilogy(1:Ncdt,Pn,'ro');
% grid on
% Ncd=Ncdt;
% ASA IRS
phi_AoD_Iw=zeros(Ncd,Nray);
theta_ZoD_Iw=zeros(Ncd,Nray);
phi_IRSA_Iw=zeros(Ncd,Nray);
theta_IRSA_Iw=zeros(Ncd,Nray);
Pn_I=zeros(Ncd,Nray);
for n0=1:Ncd
    p0=Pn(n0);
    alpha_AoD=2-4.*rand(1,Nray);
    alpha_ZoD=2-4.*rand(1,Nray);
    alpha_IRSA=2-4.*rand(1,Nray);
    alpha_IRSZ=2-4.*rand(1,Nray);
    tau=rand(1,Nray)*2;
    p_tmp=exp(-tau).*exp(-abs(alpha_AoD./c_ASD)*sqrt(2)).*exp(-abs(alpha_ZoD./c_ZSD)*sqrt(2))...
        .*exp(-abs(alpha_IRSA./c_A_IRS)*sqrt(2)).*exp(-abs(alpha_IRSZ./c_Z_IRS)*sqrt(2));
    p_tmp=p_tmp./sum(p_tmp);
    Pn_I(n0,:)=p_tmp.*p0;
    %%
    tmp_phi=alpha_AoD.*c_ASD+phi_AoD(n0); 
    tmp_phi=tmp_phi(randperm(Nray)); 
    phi_AoD_Iw(n0,:)=tmp_phi;
    tmp_phi=alpha_ZoD.*c_ZSD+theta_ZoD(n0); 
    tmp_phi=tmp_phi(randperm(Nray));
    theta_ZoD_Iw(n0,:)=tmp_phi;
    tmp_phi=alpha_IRSA.*c_A_IRS+phi_AoA(n0); 
    tmp_phi=tmp_phi(randperm(Nray)); 
    phi_IRSA_Iw(n0,:)=tmp_phi;
    tmp_phi=alpha_IRSZ.*c_Z_IRS+theta_ZoA(n0); 
    tmp_phi=tmp_phi(randperm(Nray));
    theta_IRSA_Iw(n0,:)=tmp_phi;
end  
%% user location
K=8;
xlen=5;
ylen=5;
x_location=xlen*rand(K,1);
y_location=ylen*rand(K,1);
Du=zeros(K,K);
for i0=1:K
    ix=x_location(i0);
    iy=y_location(i0);
    for j0=1:K
        jx=x_location(j0);
        jy=y_location(j0);
        dx=abs(ix-jx);
        dy=abs(iy-jy);
        Du(i0,j0)=sqrt(dx^2+dy^2);
    end
end
Tau_AD=lightspeed*2*10^(DS_mu+DS_std);
C_tau=exp(-Du./Tau_AD);
C_tau=C_tau^(1/2);
C_angle=exp(-Du./50);
C_angle=C_angle^(1/2);
%% BS-to-users
% Ncd=25;
Ncd=19;
tau_p=rand(Ncd,K).*2.*10^(DS_mu+DS_std);
tau_p=tau_p*C_tau;
phi_AoD=2*10^(ASD_mu+ASD_std).*(1-2.*rand(Ncd,K));
phi_AoD=phi_AoD*C_angle;
theta_ZoD=2*10^(ZSD_mu+ZSD_std).*(1-2.*rand(Ncd,K));
theta_ZoD=theta_ZoD*C_angle;
Pn_bu=zeros(Ncd,K);
for k0=1:K
    tau_k=tau_p(:,k0);
    tau_k=tau_k-min(tau_k);
    Zn=SFc.*randn(Ncd,1);
    Pn_p=exp(-tau_k./DS).*10.^(-Zn./10).*exp(-abs(phi_AoD(:,k0)./ASD)*sqrt(2)).*exp(-abs(theta_ZoD(:,k0)./ZSD)*sqrt(2));
%     Pn_p=exp(-tau_k./DS).*10.^(-Zn./10).*exp(-abs(phi_AoD(:,k0)./ASD)*sqrt(2));
    Pn=Pn_p./sum(Pn_p);
    Pn_bu(:,k0)=Pn;
%     [Pn, Ncb] = cluster_elimination(Pn,Thr);
end
Pn_d=zeros(Ncd,Nray,K);
phi_AoD_wu=zeros(Ncd,Nray,K);
theta_ZoD_wu=zeros(Ncd,Nray,K);
for k0=1:K
    for n0=1:Ncd
        p0=Pn_bu(n0,k0);
        alpha_AoD=2-4.*rand(1,Nray);
        alpha_ZoD=2-4.*rand(1,Nray);
        tau=rand(1,Nray)*2;
        p_tmp=exp(-tau).*exp(-abs(alpha_AoD./c_ASD)*sqrt(2)).*exp(-abs(alpha_ZoD./c_ZSD)*sqrt(2));
        p_tmp=p_tmp./sum(p_tmp);
        Pn_d(n0,:,k0)=p_tmp.*p0;
        %%
        tmp_phi=alpha_AoD.*c_ASD+phi_AoD(n0,k0); 
        tmp_phi=tmp_phi(randperm(Nray)); 
        phi_AoD_wu(n0,:,k0)=tmp_phi;
        tmp_phi=alpha_ZoD.*c_ZSD+theta_ZoD(n0,k0); 
        tmp_phi=tmp_phi(randperm(Nray));
        theta_ZoD_wu(n0,:,k0)=tmp_phi;
    end  
end
%% IRS-to-users
% Ncd=25;
Ncd=19;
tau_p=rand(Ncd,K).*2.*10^(DS_mu+DS_std);
tau_p=tau_p*C_tau;
phi_AoD=2*10^(IAS_mu+IAS_std).*(1-2.*rand(Ncd,K));
phi_AoD=phi_AoD*C_angle;
theta_ZoD=2*10^(IZS_mu+IZS_std).*(1-2.*rand(Ncd,K));
theta_ZoD=theta_ZoD*C_angle;
Pn_iu=zeros(Ncd,K);
for k0=1:K
    tau_k=tau_p(:,k0);
    tau_k=tau_k-min(tau_k);
    Zn=SFc.*randn(Ncd,1);
    Pn_p=exp(-tau_k./DS).*10.^(-Zn./10).*exp(-abs(phi_AoD(:,k0)./ASI)*sqrt(2)).*exp(-abs(theta_ZoD(:,k0)./ZSI)*sqrt(2));
%     Pn_p=exp(-tau_k./DS).*10.^(-Zn./10).*exp(-abs(phi_AoD(:,k0)./ASD)*sqrt(2));
    Pn=Pn_p./sum(Pn_p);
    Pn_iu(:,k0)=Pn;
%     [Pn, Ncb] = cluster_elimination(Pn,Thr);
end
Pn_Iu=zeros(Ncd,Nray,K);
phi_AoD_iu=zeros(Ncd,Nray,K);
theta_ZoD_iu=zeros(Ncd,Nray,K);
for k0=1:K
    for n0=1:Ncd
        p0=Pn_iu(n0,k0);
        alpha_AoD=2-4.*rand(1,Nray);
        alpha_ZoD=2-4.*rand(1,Nray);
        tau=rand(1,Nray)*2;
        p_tmp=exp(-tau).*exp(-abs(alpha_AoD./c_A_IRS)*sqrt(2)).*exp(-abs(alpha_ZoD./c_Z_IRS)*sqrt(2));
        p_tmp=p_tmp./sum(p_tmp);
        Pn_Iu(n0,:,k0)=p_tmp.*p0;
        %%
        tmp_phi=alpha_AoD.*c_A_IRS+phi_AoD(n0,k0); 
        tmp_phi=tmp_phi(randperm(Nray)); 
        phi_AoD_iu(n0,:,k0)=tmp_phi;
        tmp_phi=alpha_ZoD.*c_Z_IRS+theta_ZoD(n0,k0); 
        tmp_phi=tmp_phi(randperm(Nray));
        theta_ZoD_iu(n0,:,k0)=tmp_phi;
    end  
end
%%
save('system_para_snap_elaa.mat','fc','lightspeed','Pn_I','phi_AoD_Iw','theta_ZoD_Iw','phi_IRSA_Iw','theta_IRSA_Iw',...
    'K','x_location','y_location','Pn_d','phi_AoD_wu','theta_ZoD_wu',...
    'Pn_Iu','phi_AoD_iu','theta_ZoD_iu');

