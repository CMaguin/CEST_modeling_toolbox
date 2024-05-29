function CALC=calc_alt_var(P)
%Calculates variables for simulation computing and stores them in CALC
%structure. This lightens the computation & storage of variables
    
    %Equivalent Frequencies
    CALC.w_ref           = 2*pi*P.FREQ; %this is in 1e6 rad/s
    CALC.w1              = P.B1*gamma_2pi; %this is in rad/s
    
    %R1 and R2
    CALC.R1A=1/P.T1_water;
    CALC.R2A=1/P.T2_water;
    if P.MT
        CALC.R1MT=1/P.T1_MT;
        CALC.R2MT=1/P.T2_MT;
    end
    CALC.R1pools=1./P.T1;
    CALC.R2pools=1./P.T2;

    %Chemical shifts from resoannce frequencies
    CALC.da   = (P.xZspec-P.dw_water)*CALC.w_ref; %this is in rad/s (?)
    if P.MT; CALC.dMT = (P.xZspec-P.dw_MT)*CALC.w_ref; end
    d_pools=zeros(P.n_cest_pools,numel(P.xZspec));
    for j=1:P.n_cest_pools
        d_pools(j,:)=(P.xZspec-P.dw(j))*CALC.w_ref;
    end
    CALC.d_pools=d_pools;

    %Flip angle
    CALC.theta           = atan(CALC.w1./CALC.da);
    
    %Proton fraction
    CALC.fH=P.H ./ P.H_water;
    if P.MT; CALC.fH_MT=P.H_MT ./ P.H_water; end

    %Back exchange rates
    CALC.kbex=P.kex.*CALC.fH;
    if P.MT; CALC.kbex_MT=P.kex_MT.*CALC.fH_MT; end

    %Intramolecular exchanges!!
    if isfield(P, 'intramol_transfer_modeling')
        if P.intramol_transfer_modeling==1
            for i=1:P.n_cest_pools
                for j=1:P.n_cest_pools
                    if i==j
                        CALC.bkex_intramol(i,j)=0;
                    else
                        CALC.bkex_intramol(i,j)=P.kex_intramol(i,j)*CALC.fH(i)/CALC.fH(j);
                    end
                end
            end
        end
    end

end