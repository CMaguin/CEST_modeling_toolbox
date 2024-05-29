function A=assemble_full_matrix(Lw, Li, Kia, Kai,P, L_MT)


A=zeros(3*(1+P.n_cest_pools)+P.MT);

A(1:3,1:3)=Lw;

for i=1:P.n_cest_pools
    A(1:3,1:3)=A(1:3,1:3)-squeeze(Kai(i,:,:));
    A(1+3*i:3*(i+1),1+3*i:3*(i+1))=squeeze(Li(i,:,:)-Kia(i,:,:)); %diagonal elements
    A(1+3*i:3*(i+1),1:3)=+squeeze(Kai(i,:,:)); %"Row" i column 1 3x3 element
    A(1:3,1+3*i:3*(i+1))=squeeze(Kia(i,:,:)); %Row 1 column i 3x3 element 
end

if P.MT
    A(3,3)=A(3,3)-P.CALC.kbex_MT;
    A(end, end)=L_MT;
    A(end, 3)=P.CALC.kbex_MT;
    A(3,end)=+P.kex_MT;
end

%Intramolecular exchanges!!
do_intra=0;
if isfield(P, 'intramol_transfer_modeling')
    if P.intramol_transfer_modeling==1
        do_intra=1;
    end
end

if do_intra

    for i=1:P.n_cest_pools
    %first line (for bulk water) remains unchanged
        for j=1:P.n_cest_pools
            if j==i
                continue;
            else
                Kij=P.kex_intramol(i,j)*eye(3,3);
                Kji=P.CALC.bkex_intramol(i,j)*eye(3,3);

                %add intramolecular back-exchanges on diagonal elements
                A(1+3*i:3*(i+1),1+3*i:3*(i+1))=A(1+3*i:3*(i+1),1+3*i:3*(i+1)) -Kij; 

                %add intramolecular exchanges with pool j
                A(1+3*i:3*(i+1),1+3*j:3*(j+1))=+Kji; %Row i column j 3x3 element
            end
        end
    end
end
    
end