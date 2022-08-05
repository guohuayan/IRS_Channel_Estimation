function [dist] = distance(p1,p2)
    dist=sqrt(sum(abs(p1-p2).^2));
end

