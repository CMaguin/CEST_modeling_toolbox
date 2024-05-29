function Sol=BM_solution(P, w1, tsat,M0)

    if nargin<4
        M0=init_magnetisation_vector(P);
    end
    if nargin<3
        w1=P.CALC.w1;
        tsat=P.tp;
    end

    %Constructing the basic parts of the matrix problem
    [Lw,Li, Kia, Kai,C]=BM_matrix_parts(P,w1);
    
    Sol=zeros(3*(1+P.n_cest_pools)+P.MT,numel(P.xZspec));
    %For each offset, calculate solution
    for k=1:numel(P.xZspec)

        [Lw,Li,L_MT]=adjust_offset(Lw,Li, P, P.xZspec(k),w1);

        A=assemble_full_matrix(Lw, Li, Kia, Kai,P,L_MT);
        
        %Equation is M'=A*M+C
        %Solution is M=expm(A*tsat)*(M0+A-1*C) - A-1*C
        
        %Calculating A-1*C
        AinvC=A\C;
        
        expo_term=expm(A*tsat);
        
        Sol(:,k)=real( expo_term*(M0(:,k)+AinvC) - AinvC );
        
    end
%     
%     figure(1259)
%     plot(P.xZspec, squeeze(Sol(3,:))); hold on
    