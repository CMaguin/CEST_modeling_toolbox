function PTR=calc_PTR(f, kex, T1_water, T1_m, T2_water, T2_m, w1, tsat)
    
    R1=1/T1_water;
    p=(1+T2_water*f*kex+T2_m*kex)/(T2_m+T2_m*T2_water*f*kex);
    q=(1+T1_water*f*kex+T1_m*kex)/(T1_m+T1_water*T1_m*f*kex);

    alpha=w1^2/(w1^2+p*q);


    PTR=alpha*f*kex/(R1+f*kex)  *(1- exp(-(R1+f*kex)*tsat));



end