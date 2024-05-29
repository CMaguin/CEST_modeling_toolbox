function M=analytic_simulation(P, skipcorr)

    if nargin<2
        skipcorr=0;
    end
    
    %make sure that xZspec is a row vector
    if iscolumn(P.xZspec)
        P.xZspec=P.xZspec';
    end

    %assign "pseudo code" variables to array values
    P=convert_pc2p_var(P);
    
    if skipcorr==0
        %B1 correction? Careful to only do this once in the code
        if isfield(P,'B1c')
        P.B1=P.B1*P.B1c;
        end

        %Adjust resonance frequencies if water is not correctly centered
        %CAREFUL NOT MORE THAN ONCE IN THE CODE
        P.dw=P.dw+P.dw_water;
    end

    %Apply constrained relations, if there are some
    if isfield(P,'relations')
    for i=1:numel(P.relations)
        eval(P.relations{i});
    end
    end
    
    %Calculate alternative variables for computation
    P.CALC=calc_alt_var(P);

    %calculate R1p 
    S = calc_R1rho(P);
    M.S=S;       
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
    Zi = (P.Zi-1).*exp(-(P.CALC.R1A*(P.Trec-P.td))) +1;

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
        M.zspec=Z_cw';
    
    else %solution for pulsed sequence : see Santyr 1994 and Zaiss/Bachert 2012
        Zss         = 1 - ((1-exp(R1rho(:)*P.tp)) .* (1-cos(theta(:)) .* P.CALC.R1A ./ (-R1rho(:)))) ./ (1-exp(R1rho(:).*P.tp-P.CALC.R1A*P.td));
        Z_santyr    = (Zi-Zss).*exp((R1rho(:)*P.tp-P.CALC.R1A*P.td)*P.n) +Zss;
        M.zspec = Z_santyr;
    end
    
    if P.play_readout
        %add relaxation before readout
        if isfield(P,'readout_delay')
            M.zspec=(M.zspec-1).*exp(-(P.readout_delay/P.T1_water)) +1;
        end
        if ~isfield(P,'TEeff')
            P.TEeff=0.03;
        end
        M.zspec=(M.zspec-1).*exp(-(P.TEeff/P.T1_water)) +1;
    end
    

    if isfield(P,'Offsets_order')
    if strcmp(P.Offsets_order,'successive')
        Zi_calc=zeros(numel(P.xZspec),1);
        if isfield(P,'DS_ppm')
            P.Zi=calculate_analytic_sol(P, S, 1.0, P.DS_ppm);
            %disp(P.Zi)
        end
        for k=1:numel(P.xZspec)
            if k==1
                Zi_tmp=P.Zi;
            end
            Zi_tmp=calculate_analytic_sol(P, S, Zi_tmp, P.xZspec(k));
            Zi_calc(k)=Zi_tmp;
        end
        M.zspec=Zi_calc;
        M.Zi_calc=(1-Zi_calc).*(1-exp(-(P.CALC.R1A*P.Trec))) +Zi_calc;
    end
    end
    
    if isfield(P,'normalized')
        M0=calculate_analytic_sol(P, S, P.Zi, P.normalized);
        M.M0=M0;
        M.zspec=M.zspec./M0;
    end
    
    
    if isfield(P,'WaterResCorr')
        if P.WaterResCorr==1
            M.zspec_uncorr=M.zspec;
            if isfield(P,'WaterResCorr_limtointerp')
                corr_val=find(abs(P.xZspec)<P.WaterResCorr_limtointerp);
            else
                corr_val=find(abs(P.xZspec)<0.25);
            end
            
            Ptmp=P;
            if isfield(P,'WaterResCorr_rangetointerp')
                Ptmp.xZspec=P.WaterResCorr_rangeinterp;
            elseif isfield(P,'WaterResCorr_limtointerp')
                Ptmp.xZspec=[-1.5:0.1:-P.WaterResCorr_limtointerp,P.WaterResCorr_limtointerp:0.1:1.5];
            else
                Ptmp.xZspec=[-1.5:0.1:-0.2,0.2:0.1:1.5];
            end
            Ptmp.WaterResCorr=0;
            xtmp=[Ptmp.xZspec,P.dw_water];
            Ztmp=[analytic_simulation(Ptmp,1).zspec;0];
            
            for k=corr_val
                M.zspec(k)=interp1(xtmp,Ztmp,P.xZspec(k));
            end
        end
    end
    
    if isfield(P,'M0corr')
        M.zspec=M.zspec/P.M0corr;
    end