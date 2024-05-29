function [A,C,Lw,Li]=BM_matrix_with_blanks(P)

do_intra=0;
if isfield(P, 'intramol_transfer_modeling')
    if P.intramol_transfer_modeling==1
        do_intra=1;
    end
end

if do_intra
    
else

A=zeros(3*(1+P.n_cest_pools+P.MT));
C=zeros(3*(1+P.n_cest_pools+P.MT),1);

%Water pool
Lw=[ -P.CALC.R2A, NaN , 0;
            NaN, -P.CALC.R2A, P.CALC.w1;
            0,   -P.CALC.w1,    -P.CALC.R1A];
        
A(1:3,1:3)=Lw;
C(1)=P.CALC.R1A;

for i=1:P.n_cest_pools

    Li(i,:,:)=[ -P.CALC.R2pools(i), NaN , 0;
            NaN, -P.CALC.R2pools(i), P.CALC.w1;
            0,   -P.CALC.w1,    -P.CALC.R1pools(i)]; 

    Kai=P.CALC.kbex(i)*eye(3,3);
    Kia=P.kex(i)*eye(3,3);


    A(1:3,1:3)=A(1:3,1:3)-Kai;
    A(1+3*i:4*i,1+3*i:4*i)=Li(i,:,:)-Kia; %diagonal elements
    A(1+3*i:4*i,1)=+Kai; %"Row" i 3x3 element
    A(1,1+3*i:4*i)=Kia; %Row 1 column i 3x3 element
    
    C(1+i*3)=P.CALC.R1pools(i)*P.CALC.fH(i);
    
    
end

P.C=P.C.*P.Zi;

end
