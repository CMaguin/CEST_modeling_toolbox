function [MTR,new_x]=calc_MTRasym_with_interpolation(Z,xZspec,npoints)
    if nargin<3
        npoints=500;
    end
    densex=[linspace(min(xZspec),-0.01,floor(npoints/2)), 0, linspace(0.01,max(xZspec),floor(npoints/2)) ];
    new_x=densex;
    if isvector(Z)
        Zint=spline(xZspec, Z, densex);
        MTR=Zint(end:-1:1)-Zint;
    else
        [NbrZ, ~]=size(Z);
        MTR=zeros(NbrZ, numel(densex));
        for i=1:NbrZ
            Zint=spline(xZspec, Z(i,:), densex);
            MTR(i,:)=Zint(end:-1:1)-Zint;
        end
        
    end
    
end