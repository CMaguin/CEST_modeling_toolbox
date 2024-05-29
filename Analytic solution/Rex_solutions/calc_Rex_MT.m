function [Rex_MT, R1obs]=calc_Rex_MT(P, Reff,w1)

da      = P.CALC.da;
dc      = P.CALC.dMT;
kCA     = P.kex_MT;
kAC =   P.CALC.kbex_MT;  
r1a     = P.CALC.R1A+Reff;
r2a     = P.CALC.R2A+Reff;
r1c     = P.CALC.R1MT+Reff;

%Calculate MT lineshape
% if isfield(P, 'MT_cutoff_ppm')
%      P.MT_cutoff=P.CALC.w_ref*P.MT_cutoff_ppm; %cutoff around +-1ppm by default
% end
if ~isfield(P, 'MT_cutoff')
        P.MT_cutoff=P.CALC.w_ref; %cutoff around +-1ppm by default
end
if strcmp(P.MT_lineshape,'Flexible')
    rfmt=w1.^2.*pi.*P.MT_flex_fct(dc); %+RF_MT(P.T2_MT,w1,dc,'SuperLorentzian',P.MT_cutoff)
else
    rfmt    = RF_MT(P.T2_MT,w1,dc,P.MT_lineshape,P.MT_cutoff);
end
 
Rex_MT =    -(((da.^2 + r2a.^2).*(kCA.*r1a + (kAC + r1a).*(r1c + rfmt)) + ...
            r2a.*(kCA + r1c + rfmt).*w1.^2)./(da.^2.*(kAC + kCA + r1a + r1c + rfmt) + ...
            r2a.*(kCA.*(2.*r1a + r2a) + r2a.*(r1c + rfmt) + ...
            kAC.*(2.*r1c + r2a + 2.*rfmt) + ...
            r1a.*(2.*r1c + r2a + 2.*rfmt)) + (kCA + r1c + r2a + rfmt).*w1.^2));

R1obs =    0.5*( kAC + kCA + P.CALC.R1A+ P.CALC.R1MT - sqrt(( kAC + kCA + P.CALC.R1A + P.CALC.R1MT )^2 - 4*( kCA*P.CALC.R1A + kAC*P.CALC.R1MT + P.CALC.R1A*P.CALC.R1MT )));

end



        

        
        
 

        
