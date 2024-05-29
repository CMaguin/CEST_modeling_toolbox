function P=update_P(P)

    %assign "pseudo code" variables to array values
    P=convert_pc2p_var(P);


    
    %B1 correction?
%     if isfield(P,'B1c')
%     P.B1=P.B1*P.B1c;
%     end

    %Apply constrained relations, if there are some
    if isfield(P,'relations')
    for i=1:numel(P.relations)
        eval(P.relations{i});
    end
    end
    
    %Calculate alternative variables for computation
    P.CALC=calc_alt_var(P);

    %MT interpolation range if the case - conversion ppm to rad/s
    if isfield(P, 'MT_cutoff_ppm')
     P.MT_cutoff=P.CALC.w_ref*P.MT_cutoff_ppm; 
    end
end