function OPTIM=add_global_var(OPTIM, varname, start, lb, ub)


OPTIM.global_vars{numel(OPTIM.global_vars)+1}=varname;
OPTIM.global_start=[OPTIM.global_start start];
OPTIM.global_lb=[OPTIM.global_lb lb];
OPTIM.global_ub=[OPTIM.global_ub ub];


end