function OPTIM=add_indiv_var(OPTIM, varname, start_rel, lb_rel, up_rel)

OPTIM.indiv_vars{numel(OPTIM.indiv_vars)+1}=varname;
OPTIM.indiv_start{numel(OPTIM.indiv_start)+1}=start_rel;
OPTIM.indiv_lb{numel(OPTIM.indiv_lb)+1}=lb_rel;
OPTIM.indiv_ub{numel(OPTIM.indiv_ub)+1}=up_rel;


end