function [MTR,new_x]=calc_MTRasym_from_ref(xZspec, Z, xZspecref, Zref)
    
    [x0, ~]=find_nearest(xZspec,0);
    if isvector(Z)
        adaptedZref=spline(xZspecref, Zref, xZspec);
        noffsets=numel(Z);
        t=min(x0-1, noffsets-x0);
        Z=Z(x0-t:x0+t);
        adaptedZref=adaptedZref(x0-t:x0+t);
        new_x=xZspec(x0-t:x0+t);
        MTR=Z(end:-1:1)-adaptedZref;
    else
        original_size=size(Z);
        noffsets=original_size(end);
        sref=size(Zref);
        noffsetsref=sref(end);
        Zlin=reshape(Z, [],noffsets);
        Zlinref=reshape(Zref, [],noffsetsref);
        [NbrZ, noffsets]=size(Zlin);
        MTR=zeros(NbrZ, noffsets);
        t=min(x0-1, noffsets-x0);
        new_x=xZspec(x0-t:x0+t);
        for i=1:NbrZ
            adaptedZref=spline(xZspecref, squeeze(Zlinref(i,:)), xZspec);
            adaptedZref=adaptedZref(x0-t:x0+t)';
            Zlin(i,:)=Zlin(i,x0-t:x0+t);
            MTR(i,:)=Zlin(i,end:-1:1)-adaptedZref;
        end
        MTR=reshape(MTR, original_size);
    end
