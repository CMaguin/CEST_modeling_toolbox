function Z=calculate_analytic_sol(P, S, Zi, offset)

if nargin<3
    Zi=P.Zi;
end

if nargin>3
    P.xZspec=offset;
    P.CALC=calc_alt_var(P);
    S = calc_R1rho(P);
end

    R1rho   = S.R1rho;        
    theta   = S.theta;

%depending on the saturation type, assign Pz and Pzeff
    if strcmp(P.shape,'SL')
        Pzeff=1;
        Pz=1;  
    else
        Pzeff=cos(theta);
        Pz=cos(theta);  
    end
    


    %Correct initial magnetisation
%     Zi = (Zi-1).*exp(-(P.CALC.R1A*(P.Trec-P.td))) +1;
    Zi = (1-Zi).*(1-exp(-(P.CALC.R1A*P.Trec))) +Zi;

    %calculate solution
    if P.pulsed == 0 %continuouswave irradiation solution
        if P.MT
%             R1obs=S.R1obs;
            Z_cw = (Zi+cos(theta).*Pzeff.*S.R1obs./R1rho).*exp((R1rho*P.tp)) -cos(theta).*Pzeff.*S.R1obs./R1rho;
        else
            Z_cw = (Zi+cos(theta).*Pzeff.*P.CALC.R1A./R1rho).*exp((R1rho*P.tp)) -cos(theta).*Pzeff.*P.CALC.R1A./R1rho;
%             R1obs=P.CALC.R1A;
        end
%         Zss=cos(theta).*Pz.*R1obs./R1rho; %steady state solution
%         Z_cw= (Pz.*Pzeff.*Zi -  Zss).*exp(-(R1rho.*P.tp)) +Zss;
        Z=Z_cw';
    
    else %solution for pulsed sequence : see Santyr 1994 and Zaiss/Bachert 2012
        Zss         = 1 - ((1-exp(R1rho(:)*P.tp)) .* (1-cos(theta(:)) .* P.CALC.R1A ./ (-R1rho(:)))) ./ (1-exp(R1rho(:).*P.tp-P.CALC.R1A*P.td));
        Z_santyr    = (Zi-Zss).*exp((R1rho(:)*P.tp-P.CALC.R1A*P.td)*P.n) +Zss;
        Z = Z_santyr;
    end
    