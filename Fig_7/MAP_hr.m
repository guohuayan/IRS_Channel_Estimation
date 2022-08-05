function hat_hr = MAP_hr(hat_hr,tp_G,hat_Hd,Y1,Y2_w,Y3,bar_x,K,L2,L3,delta,iCu_w,phi2,phi3,theta0)
    for k0=1:K
        d2=(tp_G*diag(theta0))'*(Y1(:,k0)-hat_Hd(:,k0));
        for l0=1:L2
            theta=phi2(:,l0);
            d2=d2+(tp_G*diag(theta))'*(Y2_w(:,k0,l0)-hat_Hd(:,k0));
        end
        for j0=1:L3
            theta=phi3(:,j0);
            tmp1=0;
            for ii0=1:K
                tmp1=tmp1+bar_x(ii0)*hat_Hd(:,ii0);
            end
            tmp2=0;
            for ii0=1:K
                if ii0~=k0
                    tmp2=tmp2+bar_x(ii0)*tp_G*diag(theta)*hat_hr(:,ii0);
                end
            end
            d2=d2+(bar_x(k0)*tp_G*diag(theta))'*(Y3(:,j0)-tmp1-tmp2);
        end
        D2=delta^2.*iCu_w(:,:,k0)+(tp_G*diag(theta0))'*tp_G*diag(theta0);
%         D2=(tp_G*diag(theta0))'*tp_G*diag(theta0);
        for l0=1:L2
            theta=phi2(:,l0);
            D2=D2+(tp_G*diag(theta))'*tp_G*diag(theta);
        end
        for l0=1:L3
            theta=phi3(:,l0);
            D2=D2+(tp_G*diag(theta))'*tp_G*diag(theta)*abs(bar_x(k0))^2;
        end
        hat_hr(:,k0)=D2\d2;
    end
end

