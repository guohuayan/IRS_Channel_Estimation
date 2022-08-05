function s = singular(H)
    [~,s,~]=svd(H);
    s=diag(s).';
end

