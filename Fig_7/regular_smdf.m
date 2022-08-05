function A = regular_smdf(A)
    A=(A+A')/2;
    d=diag(A);
    A=A-d+real(d);
end

