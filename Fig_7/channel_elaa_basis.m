close all
clear all

load('system_para_snap_elaa.mat','fc','lightspeed','Pn_I','phi_AoD_Iw','theta_ZoD_Iw','phi_IRSA_Iw','theta_IRSA_Iw',...
    'K','x_location','y_location','Pn_d','phi_AoD_wu','theta_ZoD_wu',...
    'Pn_Iu','phi_AoD_iu','theta_ZoD_iu');
%%
rng(123);
%%
BS_loc=[0,0,3];
irs_locate=10;
xdis=10;
IRS_loc=[0,irs_locate,3];
user_loc=zeros(K,3);
for k0=1:K
    x=x_location(k0)+xdis;y=y_location(k0)+irs_locate;z=1.5;
    user_loc(k0,:)=[x,y,z];
end 
%% pathloss
Loss_BI=distance(BS_loc,IRS_loc);
Loss_BI=pathloss(Loss_BI,fc);
Loss_IU=zeros(1,K);
for k0=1:K
    p_tmp=user_loc(k0,:);
    dist=distance(p_tmp,IRS_loc);
    Loss_IU(k0)=pathloss(dist,fc);
end
Loss_BU=zeros(1,K);
for k0=1:K
    p_tmp=user_loc(k0,:);
    dist=distance(p_tmp,BS_loc);
%     Loss_BU(k0)=pathloss(dist,fc);
    Loss_BU(k0)=pathloss_nlos(dist,fc);
end
Loss_IU=Loss_BI+Loss_IU;
%% angle bias
phi_AoD_Iw=phi_AoD_Iw+90;phi_IRSA_Iw=phi_IRSA_Iw-90;
for k0=1:K
    p=user_loc(k0,:)-IRS_loc;
    [azimuth,elevation,~] = cart2sph(p(1),p(2),p(3));
    azimuth=azimuth*180/pi;
    elevation=elevation*180/pi;
    phi_AoD_iu(:,:,k0)=phi_AoD_iu(:,:,k0)+azimuth;
    theta_ZoD_iu(:,:,k0)=theta_ZoD_iu(:,:,k0)+elevation;
    p=user_loc(k0,:)-BS_loc;
    [azimuth,elevation,~] = cart2sph(p(1),p(2),p(3));
    azimuth=azimuth*180/pi;
    elevation=elevation*180/pi;
    phi_AoD_wu(:,:,k0)=phi_AoD_wu(:,:,k0)+azimuth;
    theta_ZoD_wu(:,:,k0)=theta_ZoD_wu(:,:,k0)+elevation;
end
%%
lambda=lightspeed/fc/1e9;
M1=4;M2=2;
M=M1*M2;
N1=8;N2=4;
N=N1*N2;
%% BS-IRS
% r0=1e1; %% reference distance
[Ncd,Nray]=size(phi_AoD_Iw);
Cg=zeros(M*N,M*N);
for c0=1:Ncd
%     power=Pn_I(c0,:);
%     power=power/Nray;
    for n0=1:Nray
        power=Pn_I(c0,n0);
        phi_AoD=phi_AoD_Iw(c0,n0);
        theta_AoD=theta_ZoD_Iw(c0,n0);
        alpha_B=upa_vec(M1,M2,phi_AoD,theta_AoD);
%         angle(alpha_B1)
%         alpha_B = spherical_vec(BS_array,phi_AoD,theta_AoD,r0,lambda);
%         angle(alpha_B./alpha_B(1))
%         figure
%         plot(1:8,abs(fft(alpha_B1)),'ro',1:8,abs(fft(alpha_B)),'b+')
        phi_AoA=phi_IRSA_Iw(c0,n0);
        theta_ZoA=theta_IRSA_Iw(c0,n0);
        alpha_R=upa_vec(N1,N2,phi_AoA,theta_ZoA);
%         alpha_R = spherical_vec(IRS_array,phi_AoA,theta_ZoA,r0,lambda); 
        Cg=Cg+power*kron(conj(alpha_B)*alpha_B.',alpha_R*alpha_R');
    end
end
% cond(Cg)
% [~,s,~]=svd(Cg);
% % diag(s).'
% figure
% semilogy(1:256,s,'b.');
% grid on
% TCg=zeros(N,N,M);
% for m0=1:M
%     row=(1:N)+(m0-1)*N;
%     TCg(:,:,m0)=Cg(row,row);
%     cond(TCg(:,:,m0))
% end
%% BS-User
% r0=1e5; %% reference distance
[Ncd,Nray,K]=size(phi_AoD_wu);
Cd_w=zeros(M,M,K);
for c0=1:Ncd
    for n0=1:Nray
        for k0=1:K
            power=Pn_d(c0,n0,k0);
%             power=power/Nray;
            phi_AoD=phi_AoD_wu(c0,n0,k0);
            theta_AoD=theta_ZoD_wu(c0,n0,k0);
%             alpha_B = spherical_vec(BS_array,phi_AoD,theta_AoD,r0,lambda);
            alpha_B=upa_vec(M1,M2,phi_AoD,theta_AoD);
            C_tmp=power.*alpha_B*alpha_B';
            Cd_w(:,:,k0)=Cd_w(:,:,k0)+C_tmp;
        end
    end
end
% for k0=1:K
%     cond(Cd_w(:,:,k0))
% %     rank(Cd_w(:,:,k0))
%     if rank(Cd_w(:,:,k0))~=M
%         fprintf('err rank Cd\n');
%     end
%     test=inv((Cd_w(:,:,k0)));
% end
%% IRS-user:
[Ncd,Nray,K]=size(phi_AoD_iu);
Cu_w=zeros(N,N,K);
for c0=1:Ncd
    for n0=1:Nray
        for k0=1:K
            power=Pn_Iu(c0,n0,k0);
%             power=power/Nray;
            phi_AoD=phi_AoD_iu(c0,n0,k0);
            theta_ZoD=theta_ZoD_iu(c0,n0,k0);
%             alpha_R = spherical_vec(IRS_array,phi_AoD,theta_ZoD,r0,lambda);
            alpha_R=upa_vec(N1,N2,phi_AoD,theta_ZoD);
            C_tmp=power.*alpha_R*alpha_R';
            Cu_w(:,:,k0)=Cu_w(:,:,k0)+C_tmp;
        end
    end
end
% for k0=1:K
%     cond(Cu_w(:,:,k0))
%     test=inv(Cu_w(:,:,k0));
%     if rank(Cu_w(:,:,k0))~=N
%         fprintf('err rank Cu\n');
%     end
% end
%%
Cg=regular_smdf(Cg);
for k0=1:K
    Cu_w(:,:,k0)=regular_smdf(Cu_w(:,:,k0));
    Cd_w(:,:,k0)=regular_smdf(Cd_w(:,:,k0));
end
save('elaa_covariance.mat','M','N','K','Cg','Cu_w','Cd_w','Loss_IU','Loss_BU');