function hat_G = MAP_G(hat_hr,hat_G,hat_Hd,Y1,Y2_w,Y3,bar_x,K,L2,L3,delta,iG_w,phi2,phi3,theta0,M)
    for m0=1:M
        d2=0;
%         for ii0=1:M
%             if ii0~=m0
%                 d2=d2+delta^2/2*(iG_w(:,:,m0,ii0)+iG_w(:,:,ii0,m0)')*hat_G(:,ii0);
%             end
%         end
%         d2=d2./M;
        for k0=1:K
            d2=d2+(hat_hr(:,k0).'*diag(theta0))'*(Y1(m0,k0)-hat_Hd(m0,k0));
            for l0=1:L2
                theta=phi2(:,l0);
                d2=d2+(hat_hr(:,k0).'*diag(theta))'*(Y2_w(m0,k0,l0)-hat_Hd(m0,k0));
            end
        end
        tmp1=0;
        tmp2=0;
        for ii0=1:K
            tmp1=tmp1+bar_x(ii0)*hat_hr(:,ii0).';
            tmp2=tmp2+bar_x(ii0)*hat_Hd(m0,ii0);
        end
        for j0=1:L3
            theta=phi3(:,j0);
            d2=d2+(tmp1*diag(theta))'*(Y3(m0,j0)-tmp2);
        end
        D1=delta^2*iG_w(:,:,m0,m0);
%         D1=0;
        for k0=1:K
            D1=D1+(hat_hr(:,k0).'*diag(theta0))'*hat_hr(:,k0).'*diag(theta0);
            for l0=1:L2
                theta=phi2(:,l0);
                D1=D1+(hat_hr(:,k0).'*diag(theta))'*hat_hr(:,k0).'*diag(theta);
            end
        end
        for l0=1:L3
            theta=phi3(:,l0);
            D1=D1+(tmp1*diag(theta))'*tmp1*diag(theta);
        end
        hat_G(:,m0)=D1\d2;
    end
end

