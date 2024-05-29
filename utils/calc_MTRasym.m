function [MTR,new_x]=calc_MTRasym(Z,xZspec)
    [x0, ~]=find_nearest(xZspec,0);
    if isvector(Z)
        noffsets=numel(Z);
        t=min(x0-1, noffsets-x0);
        Z=Z(x0-t:x0+t);
        new_x=xZspec(x0-t:x0+t);
        MTR=Z(end:-1:1)-Z;
    else
        original_size=size(Z);
        noffsets=original_size(end);
        Zlin=reshape(Z, [],noffsets);
        [NbrZ, noffsets]=size(Zlin);
        MTR=zeros(NbrZ, noffsets);
        t=min(x0-1, noffsets-x0);
        new_x=xZspec(x0-t:x0+t);
        for i=1:NbrZ
            Zlin(i,:)=Zlin(i,x0-t:x0+t);
            MTR(i,:)=Zlin(i,end:-1:1)-Zlin(i,:);
        end
        MTR=reshape(MTR, original_size);
    end
    
end