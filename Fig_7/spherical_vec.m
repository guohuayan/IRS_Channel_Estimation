function alpha_B = spherical_vec(BS_array,phi_AoD,theta_AoD,r0,lambda)
    phi_AoD=phi_AoD./180.*pi;
    theta_AoD=theta_AoD./180.*pi;
    [M]=size(BS_array,1);
    dist_w=zeros(M,1);
    point=sph_to_cart(phi_AoD,theta_AoD,r0);
    for m0=1:M
        point_a=BS_array(m0,:);
        dist_w(m0)=distance_cart(point_a,point);
    end
    alpha_B=exp(-1j.*2*pi./lambda.*dist_w)./dist_w;
    alpha_B=alpha_B./sqrt(mean(abs(alpha_B).^2));
end

function point=sph_to_cart(phi_AoD,theta_AoD,r0)
    point=zeros(1,3);
    point(1)=r0.*cos(phi_AoD).*sin(theta_AoD);
    point(2)=r0.*sin(phi_AoD).*sin(theta_AoD);
    point(3)=r0.*cos(theta_AoD);
end
function dist=distance_cart(p1,p2)
    dist=sqrt(sum(abs(p1-p2).^2));
end

