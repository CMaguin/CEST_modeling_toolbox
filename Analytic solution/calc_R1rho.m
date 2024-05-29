function S =calc_R1rho(P)

% if isrow(P.xZspec)
%     P.xZspec=P.xZspec';
% end

noffsets=numel(P.xZspec);
w_ref           = P.CALC.w_ref;
w1              = P.CALC.w1;
da              = P.CALC.da;
d_pools         = P.CALC.d_pools;
theta           = P.CALC.theta;
S.theta         = theta;
S.Reff_sincos   = -P.CALC.R1A*cos(theta).^2 -(P.CALC.R2A)*sin(theta).^2;

if strcmp(P.shape,'block')
    switch P.Rex_sol    
            case 'Hyper' 
                S.Rex_hyper = calc_Rex_hyper(P);
                if P.MT
                    [S.Rex_MT, S.R1obs] = calc_Rex_MT(P,S.Reff_sincos, P.CALC.w1);
                    S.R1rho = S.Reff_sincos + S.Rex_hyper./(1+P.CALC.fH_MT) + S.Rex_MT;
                else
                    S.R1rho = S.Reff_sincos + S.Rex_hyper;
                end
    
            case 'Lorentz'
                S.Rex_lorentz = calc_Rex_lorentz(P);
                if P.MT
                    [S.Rex_MT, S.R1obs] = calc_Rex_MT(P,S.Reff_sincos, P.CALC.w1);
                    S.R1rho = S.Reff_sincos + S.Rex_lorentz./(1+P.CALC.fH_MT) + S.Rex_MT;
                else
                    S.R1rho = S.Reff_sincos + S.Rex_lorentz;
                end
    
             case 'minilorentz'
                S.Rex_minilorentz = calc_Rex_minilorentz(P);
                if P.MT
                    [S.Rex_MT, S.R1obs] = calc_Rex_MT(P,S.Reff_sincos, P.CALC.w1);
                    S.R1rho = S.Reff_sincos + S.Rex_minilorentz./(1+P.CALC.fH_MT) + S.Rex_MT;
                else
                    S.R1rho = S.Reff_sincos + S.Rex_minilorentz;
                end
            otherwise
                error('Unrecognised Rex solution type')
    end


else %other shapes of solution

    integration_seg=200;
    [B1]=RF_pulse(P,integration_seg); %calculate B1 pulse depending on the pulse shape given
    
    %Preallocating for speed
    Reff_mean   = ones(1,noffsets);
    Rex_mean=ones(P.n_cest_pools,noffsets);
    
    
    %Calculate integrated Reff and Rex 
    for i=1:noffsets
        Reff_mean(i) = trapz(Reffw1(da(i),gamma_2pi*B1,P.CALC.R2A,P.CALC.R1A))/integration_seg ;
        for j=1:P.n_cest_pools
            Rex_mean(j,i)=trapz(Rex_Hyper_onepool(da(i),gamma_2pi*B1,d_pools(j,i),P.kex(j),P.CALC.kbex(j),P.CALC.R2pools(j)))/integration_seg ;
        end
    end

    %Full Rex (sum over all pools)
    S.Rex_full=sum(Rex_mean,1);
    S.R1rho=Reff_mean+S.Rex_full;

    %S.alpha_f_mean=-sum(Rex_mean./P.kex,1);

    %MT pool
    if P.MT==1
        S.R1obs     = 0.5*( P.kex_MT + P.CALC.kbex_MT + P.CALC.R1A + P.CALC.R1MT - sqrt(( P.CALC.kbex_MT + P.kex_MT + P.CALC.R1A + P.CALC.R1MT )^2 - 4*( P.kex_MT*P.CALC.R1A + P.CALC.kbex_MT*P.CALC.R1MT + P.CALC.R1A*P.CALC.R1MT )));
        for j=1:numel(B1)
            Rex_MT_full_B1(j,:)=calc_Rex_MT(P,Reffw1(da,gamma_2pi*B1(j),P.CALC.R2A,P.CALC.R1A), gamma_2pi*B1(j));
        end
        Rex_mean_MT=trapz(Rex_MT_full_B1,1)/integration_seg;
        S.Rex_mean_MT   = Rex_mean_MT;
        S.R1rho      = Reff_mean+S.Rex_full./(1+P.CALC.fH_MT)+ S.Rex_mean_MT;
        %S.alpha_f_mean = S.alpha_f_mean - Rex_mean_MT/P.kCA;
    end

end