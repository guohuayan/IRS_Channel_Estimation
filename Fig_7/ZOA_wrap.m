function [tmp_theta] = ZOA_wrap(tmp_theta)
    while(1)
        tmp_theta=tmp_theta+(tmp_theta<0)*360-(tmp_theta>360)*360;
        len=sum((tmp_theta<0))+sum((tmp_theta>360));
        if len==0
            break
        end
    end
    tmp_g1=tmp_theta<180;
    tmp_g2=1-tmp_g1;
    tmp_theta=tmp_theta.*tmp_g1+(360-tmp_theta).*tmp_g2;
end

