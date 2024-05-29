function [Lw,Li,L_MT]=adjust_offset(Lw,Li, P, drf, w1)

%Water part
wa=-(drf-P.dw_water)*P.CALC.w_ref;
Lw(1,2)=wa;
Lw(2,1)=(-wa);

for i=1:P.n_cest_pools
    wi=-(drf-P.dw(i))*P.CALC.w_ref;
    Li(i,1,2)=wi;
    Li(i,2,1)=(-wi);
end

if P.MT
    w_MT=(drf-P.dw_MT)*P.CALC.w_ref;
    L_MT=-P.CALC.R1MT-RF_MT(P.T2_MT,w1,w_MT,P.MT_lineshape)-P.kex_MT;
else
    L_MT=[];
end