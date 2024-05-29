function OPTIM=add_voxel_var(OPTIM, varname, start_rel, lb_rel, up_rel)

OPTIM.voxel_vars{numel(OPTIM.voxel_vars)+1}=varname;
OPTIM.voxel_start{numel(OPTIM.voxel_start)+1}=start_rel;
OPTIM.voxel_lb{numel(OPTIM.voxel_lb)+1}=lb_rel;
OPTIM.voxel_ub{numel(OPTIM.voxel_ub)+1}=up_rel;


end