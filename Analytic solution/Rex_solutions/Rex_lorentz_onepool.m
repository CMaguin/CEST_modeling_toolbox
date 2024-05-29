function Rex_Lorentz=Rex_lorentz_onepool(da,w1,di,ki,ka,r2i)
    % calculate Rex Lorentz CEST solution for pool i
    %NBM paper;    assumes (kAB<<kBA, R1B<<kBA, Reff<<R2B)

%     w1=w1*2*pi; da=da*2*pi; di=di*2*pi;
    REXMAX  = -ka.*w1^2./(da.^2+w1.^2).*((da-di).^2 +(da.^2+w1.^2).*r2i./ki + r2i.*(ki+r2i));
%     REXMAX  = -ka.*sin(theta).^2.*((da-di).^2 +(da.^2+w1.^2).*r2i./ki + r2i.*(ki+r2i));
    GAMMA   = 2*sqrt( (ki+r2i)./ki.*w1.^2 + (ki+r2i).^2);
    Rex_Lorentz = REXMAX./((GAMMA./2).^2+di.^2);
    
%     pf=-ka.*w1^2./(da.^2+w1.^2) ;
%     a1=(da-di).^2;
%     a2=(da.^2+w1.^2).*r2i./ki;
%     b=r2i.*(ki+r2i);
%     
%     gamma1=(ki+r2i)./ki.*w1.^2;
%     gamma2=(ki+r2i).^2;
end