function Rex_miniLorentz=Rex_minilorentz_onepool(da,w1,di,ki,ka,r2i)
% calculate Rex miniLorentz CEST solution for pool i
% w1=w1*2*pi;
% da=da*2*pi;
% di=di*2*pi;
    REXMAX  = -ka.*w1^2./(w1.^2+ki.*(ki+r2i));
    GAMMA   = 2*sqrt( (ki+r2i)./ki.*w1.^2 + (ki+r2i).^2);
    Rex_miniLorentz = REXMAX*(GAMMA./2).^2./((GAMMA./2).^2+di.^2);
end