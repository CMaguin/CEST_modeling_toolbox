function [Rex_MT, rfmt] = calc_Rex_MT2(da,w1,dc,P)
    
    theta   = P.CALC.theta;
    Reff    = -(P.CALC.R1A*cos(theta).^2 + P.CALC.R2A*sin(theta).^2);
    r1a     = P.CALC.R1A+Reff;
    r2a     = P.CALC.R2A+Reff;
    r1c     = P.CALC.R1MT+Reff;
    rfmt    = RF_MT(P.T2_MT,w1,dc,P.MT_lineshape);

    kCA     = P.kex_MT;
    kAC = P.CALC.kbex_MT;

    if isrow(rfmt)
        Rex_MT =    -(((da.^2 + r2a.^2).*(kCA.*r1a + (kAC + r1a).*(r1c + rfmt)) + ...
                    r2a.*(kCA + r1c + rfmt).*w1.^2)./(da.^2.*(kAC + kCA + r1a + r1c + rfmt) + ...
                    r2a.*(kCA.*(2.*r1a + r2a) + r2a.*(r1c + rfmt) + ...
                    kAC.*(2.*r1c + r2a + 2.*rfmt) + ...
                    r1a.*(2.*r1c + r2a + 2.*rfmt)) + (kCA + r1c + r2a + rfmt).*w1.^2));
    else
        [~, nint]=size(rfmt);
        for k=1:nint
            rfmt_indint=rfmt(:,k);
            rfmt_indint=rfmt_indint';
            R =    -(((da.^2 + r2a.^2).*(kCA.*r1a + (kAC + r1a).*(r1c + rfmt_indint)) + ...
            r2a.*(kCA + r1c + rfmt_indint).*w1.^2)./(da.^2.*(kAC + kCA + r1a + r1c + rfmt_indint) + ...
            r2a.*(kCA.*(2.*r1a + r2a) + r2a.*(r1c + rfmt_indint) + ...
            kAC.*(2.*r1c + r2a + 2.*rfmt_indint) + ...
            r1a.*(2.*r1c + r2a + 2.*rfmt_indint)) + (kCA + r1c + r2a + rfmt_indint).*w1.^2));
            Rex_MT(k,:)=R;
        end
    end
end


