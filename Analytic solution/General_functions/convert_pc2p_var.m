function P=convert_pc2p_var(P)

for i=1:P.n_cest_pools
   
    if isfield(P, strcat('kex_',P.pool_names{i}))
       P.kex(strcmp( P.pool_names, P.pool_names{i}))=P.(strcat('kex_',P.pool_names{i}));
    end
    
    if isfield(P, strcat('H_',P.pool_names{i}))
       P.H(strcmp( P.pool_names, P.pool_names{i}))=P.(strcat('H_',P.pool_names{i}));
    end
    
    if isfield(P, strcat('dw_',P.pool_names{i}))
       P.dw(strcmp( P.pool_names, P.pool_names{i}))=P.(strcat('dw_',P.pool_names{i}));
    end
    
    if isfield(P, strcat('T1_',P.pool_names{i}))
       P.T1(strcmp( P.pool_names, P.pool_names{i}))=P.(strcat('T1_',P.pool_names{i}));
    end
    
    if isfield(P, strcat('T2_',P.pool_names{i}))
       P.T2(strcmp( P.pool_names, P.pool_names{i}))=P.(strcat('T2_',P.pool_names{i}));
    end
    
end

