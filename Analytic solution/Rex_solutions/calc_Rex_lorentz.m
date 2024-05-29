function Rex=calc_Rex_lorentz(P)

    if P.n_cest_pools==0
        Rex=0;
    else
        Individual_Rex=zeros(P.n_cest_pools,numel(P.xZspec));
        for j=1:P.n_cest_pools
            Individual_Rex(j,:)=Rex_lorentz_onepool(P.CALC.da,P.CALC.w1,P.CALC.d_pools(j,:),P.kex(j),P.CALC.kbex(j),P.CALC.R2pools(j));
        end
        Rex=sum(Individual_Rex,1);
        
%         figure(846)
%         plot(P.xZspec, Individual_Rex); hold on
%         plot(P.xZspec, Rex, 'b')
    end
end


