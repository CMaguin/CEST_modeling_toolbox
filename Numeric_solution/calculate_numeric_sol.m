function M_sol=calculate_numeric_sol(P, Mi, offset)

if nargin<2
    Mi=P.Zi;
end

if nargin>2
    P.xZspec=offset;
    P.CALC=calc_alt_var(P);
end

    %recovery from Zi
    M_sol=BM_solution(P,0,P.Trec,Mi);
    
    %Saturation sequence
    if P.pulsed==0
        M_sol=BM_solution(P,P.CALC.w1,P.tp,M_sol);
    elseif strcmp(P.shape,'block')
        for i=1:P.n
            %saturation
            M_sol=BM_solution(P,P.CALC.w1,P.tp,M_sol);
            
            %relaxation during td
            %%%M_sol=(M_sol-1).*exp(-(P.CALC.R1A*P.td)) +1;
            M_sol=BM_solution(P,0,P.td,M_sol);
        end
        

        
    else
        %pulse decomposition
        steps=50; 
        time_interval=P.tp/steps;
        [B1]=RF_pulse(P, steps);
        
        %init magnetization vector
        M_sol=init_magnetisation_vector(P);
        for n=1:P.n %loop on each pulse
            %saturation pulse
            for i=1:steps
               M_sol=BM_solution(P,gamma_2pi*B1(i),time_interval,M_sol);
            end
            %relaxation
            M_sol=BM_solution(P,0,P.td,M_sol);
        end
    end
    
    %Readout phase = fast spin echo sequence 
    if P.play_readout
        x_index=find(mod(1:3*(1+P.n_cest_pools)+P.MT, 3)==0);
        y_index=find(mod(1:3*(1+P.n_cest_pools)+P.MT, 3)==1);
        z_index=find(mod(1:3*(1+P.n_cest_pools)+P.MT, 3)==2);
        
        if ~isfield(P,'flipangle'); P.flipangle=90; end
        if ~isfield(P,'TEeff'); P.TEeff=0.03; end
        if ~isfield(P,'spoilf'); P.spoilf=0; end %0=full spoiling; 1=no spoiling
        %if ~isfield(P,'RAREfactor'); P.RAREfactor=90; end
        %t_rf=0.005;
        %B1rf=P.flipangle/(360*gamma_*t_rf);
        %M_sol(y_index)=cosd(flipangle)*M_sol(y_index)+sind(flipangle)*M_sol(z_index);
        %M_sol(z_index)=-sind(flipangle)*M_sol(y_index)+cosd(flipangle)*M_sol(z_index);
        M_sol(x_index)=P.spoilf*M_sol(x_index); M_sol(y_index)=P.spoilf*M_sol(y_index); %spoiling
        M_sol=BM_solution(P,0,P.TEeff,M_sol); %effective echo time relaxation
    end
    