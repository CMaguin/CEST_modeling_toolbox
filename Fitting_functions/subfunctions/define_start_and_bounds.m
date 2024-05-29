function OPTIM=define_start_and_bounds(OPTIM,P)

OPTIM.start=[];
OPTIM.lowerbounds=[];
OPTIM.upperbounds=[];

%First global values
for k=1:numel(OPTIM.global_vars)
    OPTIM.start=[OPTIM.start; OPTIM.global_start(k)];
    OPTIM.lowerbounds=[OPTIM.lowerbounds; OPTIM.global_lb(k)];
    OPTIM.upperbounds=[OPTIM.upperbounds; OPTIM.global_ub(k)];
end

%Then voxel values
voxel_start=zeros(numel(OPTIM.voxel_vars)*OPTIM.nvoxels,1); voxel_lowerbounds=zeros(numel(OPTIM.voxel_vars)*OPTIM.nvoxels,1); voxel_upperbounds=zeros(numel(OPTIM.voxel_vars)*OPTIM.nvoxels,1);
indiv_start=zeros(numel(OPTIM.indiv_vars)*OPTIM.nvoxels*OPTIM.nconditions,1); indiv_lowerbounds=zeros(numel(OPTIM.indiv_vars)*OPTIM.nvoxels*OPTIM.nconditions,1);indiv_upperbounds=zeros(numel(OPTIM.indiv_vars)*OPTIM.nvoxels*OPTIM.nconditions,1);
for i=1:OPTIM.nvoxels
    for j=1:numel(OPTIM.vary_voxel)
        v=OPTIM.vary_voxel{j};
        P.(v)=OPTIM.varyval_voxel(j,i);
    end
    for k=1:numel(OPTIM.voxel_vars)
        voxel_start(sub2ind([OPTIM.nvoxels,numel(OPTIM.voxel_vars)],i,k))= eval(OPTIM.voxel_start{k});
        voxel_lowerbounds(sub2ind([OPTIM.nvoxels,numel(OPTIM.voxel_vars)],i,k))=eval(OPTIM.voxel_lb{k});
        voxel_upperbounds(sub2ind([OPTIM.nvoxels,numel(OPTIM.voxel_vars)],i,k))=eval(OPTIM.voxel_ub{k});
    end

    for icond=1:OPTIM.nconditions
        for j=1:numel(OPTIM.vary_indiv)
        v=OPTIM.vary_indiv{j};
        P.(v)=OPTIM.varyval_indiv(j,icond);
        end
        for k=1:numel(OPTIM.indiv_vars)
        indiv_start(sub2ind([OPTIM.nvoxels,OPTIM.nconditions,numel(OPTIM.indiv_vars)],i,icond,k))= eval(OPTIM.indiv_start{k});
        indiv_lowerbounds(sub2ind([OPTIM.nvoxels,OPTIM.nconditions,numel(OPTIM.indiv_vars)],i,icond,k))=eval(OPTIM.indiv_lb{k});
        indiv_upperbounds(sub2ind([OPTIM.nvoxels,OPTIM.nconditions,numel(OPTIM.indiv_vars)],i,icond,k))=eval(OPTIM.indiv_ub{k});
        end
    end
end

OPTIM.start=[OPTIM.start; voxel_start; indiv_start];
OPTIM.lowerbounds=[OPTIM.lowerbounds; voxel_lowerbounds; indiv_lowerbounds];
OPTIM.upperbounds=[OPTIM.upperbounds; voxel_upperbounds; indiv_upperbounds];

% OPTIM.start=OPTIM.start';
% OPTIM.lowerbounds=OPTIM.lowerbounds';
% OPTIM.upperbounds=OPTIM.upperbounds';

OPTIM.dep_vars=OPTIM.global_vars;



