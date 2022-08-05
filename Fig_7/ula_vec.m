function alpha = ula_vec(M,theta)
    theta=theta/180*pi;
    alpha=1:M;
    alpha=alpha-1;
    alpha=alpha.';
    alpha=1j.*pi.*sin(theta).*alpha;
    alpha=exp(-alpha);
end

