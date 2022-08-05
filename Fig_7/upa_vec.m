function alpha = upa_vec(N1,N2,phi,xi)
    phi=phi/180*pi;
    xi=xi/180*pi;
    alpha1=1:N1;
    alpha1=alpha1-1;
    alpha1=1j.*pi.*sin(xi).*cos(phi).*alpha1;
    alpha1=exp(-alpha1);
    alpha1=alpha1.';
    alpha2=1:N2;
    alpha2=alpha2-1;
    alpha2=1j.*pi.*sin(xi).*sin(phi).*alpha2;
    alpha2=exp(-alpha2);
    alpha2=alpha2.';
    alpha=kron(alpha1,alpha2);
end

