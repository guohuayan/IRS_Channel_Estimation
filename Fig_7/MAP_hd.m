function hat_Hd = MAP_hd(hat_hr,tp_G,hat_Hd,Y1,Y2_w,Y3,bar_x,K,L2,L3,delta,iCd_w,phi2,phi3,theta0,M)
    for k0=1:K
        d1=Y1(:,k0)-tp_G*diag(theta0)*hat_hr(:,k0);
        for l0=1:L2
            theta=phi2(:,l0);
            d1=d1+Y2_w(:,k0,l0)-tp_G*diag(theta)*hat_hr(:,k0);
        end
        for j0=1:L3
            theta=phi3(:,j0);
            tmp1=0;
            for ii0=1:K
                tmp1=tmp1+bar_x(ii0)*tp_G*diag(theta)*hat_hr(:,ii0);
            end
            tmp2=0;
            for ii0=1:K
                if ii0~=k0
                    tmp2=tmp2+bar_x(ii0)*hat_Hd(:,ii0);
                end
            end
            d1=d1+bar_x(k0)'*(Y3(:,j0)-tmp1-tmp2);
        end
        tmp=((L2+1+L3*abs(bar_x(k0))^2)*eye(M)+delta^2.*iCd_w(:,:,k0))\d1;
%         tmp=((L2+1+L3*abs(bar_x(k0))^2)*eye(M))\d1;
        hat_Hd(:,k0)=tmp;
    end
end

