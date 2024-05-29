function PMEX=update_pulsedqcest(PMEX,P)

%update PMEX 
PMEX.WaterPool.R1=P.CALC.R1A;
PMEX.WaterPool.R2=P.CALC.R2A;

for i=1:P.n_cest_pools
    PMEX.CESTPool(i).f=P.CALC.fH(i);
    PMEX.CESTPool(i).R1=1/P.T1(i);
    PMEX.CESTPool(i).R2=1/P.T2(i);
    PMEX.CESTPool(i).k=P.kex(i);
    PMEX.CESTPool(i).dw=P.dw(i);
end