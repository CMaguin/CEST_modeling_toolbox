function OPTIM=init_OPTIM_struc(nvoxels, nconditions)
    if nargin==0
        nvoxels=1;
        nconditions=1;
    end

    OPTIM.multivar=1;
    
    OPTIM.global_vars={};
    OPTIM.global_start=[];
    OPTIM.global_lb=[];
    OPTIM.global_ub=[];
    
    OPTIM.nvoxels=nvoxels;
    OPTIM.voxel_vars={};
    OPTIM.voxel_start={};
    OPTIM.voxel_lb={};
    OPTIM.voxel_ub={};

    OPTIM.vary_voxel={};
    OPTIM.varyval_voxel=[];
    
    OPTIM.nconditions=nconditions;
    OPTIM.indiv_vars={};
    OPTIM.indiv_start={};
    OPTIM.indiv_lb={};
    OPTIM.indiv_ub={};   

    OPTIM.vary_indiv={};
    OPTIM.varyval_indiv=[];