function M=numeric_simulation_lipoGlu(P)

    if iscolumn(P.xZspec)
        P.xZspec=P.xZspec';
    end

    %assign "pseudo code" variables to array values
    P=convert_pc2p_var(P);
    
    %B1 correction?
    if isfield(P,'B1c')
        P.B1=P.B1*P.B1c;
    end
    
    %Adjust resonance frequencies if water is not correctly centered
    P.dw=P.dw+P.dw_water;

    %Apply constrained relations, if there are some
    if isfield(P,'relations')
        for i=1:numel(P.relations)
            eval(P.relations{i});
        end
    end
    
    %Calculate alternative variables for computation
    P.CALC=calc_alt_var(P);
    
    M_sol=P.Zi*init_magnetisation_vector(P);
    
    
    
    if isfield(P,'Trec')
    %recovery from Zi
    M_sol=BM_solution_lipoGlu(P,0,P.Trec,M_sol);
    end
    
    %Saturation sequence
    if P.pulsed==0
        M_sol=BM_solution_lipoGlu(P,P.CALC.w1,P.tp,M_sol);
        M.sol=M_sol;
        M.zspec=M_sol(3,:); % take water z-magnetization as solution
    elseif strcmp(P.shape,'block')
        for i=1:P.n
            %saturation
            M_sol=BM_solution_lipoGlu(P,P.CALC.w1,P.tp,M_sol);
            
            %relaxation during td
            %%%M_sol=(M_sol-1).*exp(-(P.CALC.R1A*P.td)) +1;
            M_sol=BM_solution_lipoGlu(P,0,P.td,M_sol);
        end
        
        M.sol=M_sol;
        
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
               M_sol=BM_solution_lipoGlu(P,gamma_2pi*B1(i),time_interval,M_sol);
            end
            %relaxation
            M_sol=BM_solution_lipoGlu(P,0,P.td,M_sol);
        end
         M.sol=M_sol;    
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
        M_sol=BM_solution_lipoGlu(P,0,P.TEeff,M_sol); %effective echo time relaxation
    end
    
    if isfield(P,'Offsets_order')
    if strcmp(P.Offsets_order,'successive')
        Mi_calc=zeros(3*(1+P.n_cest_pools)+P.MT,numel(P.xZspec));
        for k=1:numel(P.xZspec)
            if k==1
                Mi_tmp=P.Zi*init_magnetisation_vector(P);
            end
            Mi_tmp=calculate_numeric_sol_lipoGlu(P, Mi_tmp, P.xZspec(k));
            Mi_calc(:,k)=Mi_tmp;
        end
        M_sol=Mi_calc;
    end
    end
    
    M.sol=M_sol;
    
    M.zspec=M_sol(3,:); % take water z-magnetization as solution
    
   if isfield(P,'normalized')
        M0=calculate_numeric_sol_lipoGlu(P, P.Zi*init_magnetisation_vector(P), P.normalized);
        M.M0=M0;
        M.zspec=M.zspec./M0(3);
    end
    
    if isrow(M.zspec)
        M.zspec=M.zspec';
    end
    
    