function Rex=calc_Rex_minilorentz(P)

    if P.n_cest_pools==0
        Rex=0;
    else
        Individual_Rex=zeros(P.n_cest_pools,1);
        for j=1:P.n_cest_pools
            Individual_Rex(j)=Rex_minilorentz_onepool(P.CALC.da,P.CALC.w1,P.CALC.d_pools(j),P.kex(j),P.CALC.kbex(j),P.CALC.R2pools(j));
        end
        Rex=sum(Individual_Rex);
    end
end


