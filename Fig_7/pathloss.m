function [path] = pathloss(d,fc)
%     path=32.4+21*log10(d)+20*log10(fc);
    path=32.4+17.3*log10(d)+20*log10(fc);
end

