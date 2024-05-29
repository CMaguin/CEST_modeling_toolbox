function [Lw,Li, Kia, Kai,C]=BM_matrix_parts(P,w1)


C=zeros(3*(1+P.n_cest_pools)+P.MT,1);

%Water pool
Lw=[ -P.CALC.R2A, NaN , 0;
            NaN, -P.CALC.R2A, w1;
            0,   -w1,    -P.CALC.R1A];
        
C(3)=P.CALC.R1A;
Li=zeros(P.n_cest_pools, 3,3);
Kai=zeros(P.n_cest_pools, 3,3);
Kia=zeros(P.n_cest_pools, 3,3);
for i=1:P.n_cest_pools

    Li(i,:,:)=[ -P.CALC.R2pools(i), NaN , 0;
                NaN, -P.CALC.R2pools(i), w1;
                0,   -w1,    -P.CALC.R1pools(i)]; 

    Kai(i,:,:)=P.CALC.kbex(i)*eye(3,3);
    Kia(i,:,:)=P.kex(i)*eye(3,3);
    
    C(3*(i+1))=P.CALC.R1pools(i)*P.CALC.fH(i);
    
    
end


if P.MT
    C(end)=P.CALC.R1MT*P.CALC.fH_MT;
end

if isfield(P, 'intramol_transfer_modeling')
    if P.intramol_transfer_modeling==1
        %need to add intramolecular transfer parts

        %warning('Didnt add inter metabolite exchange yet')
    end
end
