function y0=init_magnetisation_vector(P)

y0=zeros(3*(1+P.n_cest_pools)+P.MT,numel(P.xZspec));

%set water pool magnetisation to Zi
y0(3,:)=1.0;

for i=1:P.n_cest_pools
    y0((i+1)*3,:)=P.CALC.fH(i);
end

if P.MT
    y0(end,:)=P.CALC.fH_MT;
end

%y0=y0*P.Zi;