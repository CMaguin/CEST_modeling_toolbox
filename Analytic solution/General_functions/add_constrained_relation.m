function P=add_constrained_relation(P,rel)
    if ~isfield(P,'relations')
        P.relations={};
    end

    P.relations{numel(P.relations)+1}=rel;

end