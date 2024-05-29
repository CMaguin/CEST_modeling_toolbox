function [P, OPTIM_new]=assign_results(Results,P,OPTIM)

%Global variables
for i =1:numel(OPTIM.global_vars)
   P.(OPTIM.global_vars{i})=Results.(OPTIM.global_vars{i});
end

Nvaryvoxels=numel(OPTIM.vary_voxel);
voxel_vars=OPTIM.voxel_vars;
voxel_val=OPTIM.varyval_voxel;
for i=1:numel(Results.voxel_vars)
    if ~any(strcmp(OPTIM.vary_voxel{i}, Results.voxel_vars{i}))
        Nvaryvoxels=Nvaryvoxels+1;
        voxel_vars{numel(voxel_vars)+1}=OPTIM.voxel_vars{i};
        voxel_val(numel(voxel_vars)+1,:)=Results.voxel_val(i,:);
    else
        voxel_val(strcmp(OPTIM.vary_voxel{i}, OPTIM.voxel_vars{i}),:)=Results.voxel_val(i,:);
    end
end

Nvaryindiv=numel(OPTIM.vary_indiv);
indiv_vars=OPTIM.vary_indiv;
indiv_val=OPTIM.varyval_indiv;
for i=1:numel(Results.indiv_vars)
    if ~any(strcmp(OPTIM.vary_indiv, Results.indiv_vars{i}))
        Nvaryindiv=Nvaryindiv+1;
        indiv_vars{numel(indiv_vars)+1}=Results.indiv_vars{i};
        indiv_val(numel(indiv_vars)+1,:)=Results.indiv_val(i,:);
    else
        indiv_val(strcmp(OPTIM.vary_indiv{i}, OPTIM.indiv_vars{i}),:)=Results.indiv_val(i,:);
    end
end
        
OPTIM_new=init_OPTIM_struc(OPTIM.nvoxels, OPTIM.nconditions);

OPTIM_new.vary_voxel=voxel_vars;
OPTIM_new.varyval_voxel=voxel_val;

OPTIM_new.vary_indiv=indiv_vars;
OPTIM_new.varyval_indiv=indiv_val;

