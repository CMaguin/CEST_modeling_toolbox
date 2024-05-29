function zspec=compute_analytical_multivar(P, OPTIM, assigned_values,xdata)

zspec=zeros(OPTIM.nvoxels, OPTIM.nconditions, numel(P.xZspec));

if isfield(P, 'sim_type')
    if P.sim_type==0 
        sim=@(P) numeric_simulation(P); 
    elseif P.sim_type==-1
        sim=@(P) pulsedq_cest_simu(P);
    else
        sim=@(P) analytic_simulation(P); 
    end
else
    sim=@(P) analytic_simulation(P);
end


if nargin>3
    P.xZspec=xdata;
end

if iscolumn(P.xZspec)
    P.xZspec=P.xZspec';
end

if OPTIM.fit_type==2
    zspec=zeros(OPTIM.nvoxels, OPTIM.nconditions, 2*numel(P.xZspec));
else
    zspec=zeros(OPTIM.nvoxels, OPTIM.nconditions, numel(P.xZspec));
end


%Separate global from voxel and individual variables
ig=numel(OPTIM.global_vars);
global_values=assigned_values(1:ig);
iv=numel(OPTIM.voxel_vars)*OPTIM.nvoxels;
voxel_values=assigned_values(1+ig:ig+iv);
iin=numel(OPTIM.indiv_vars)*OPTIM.nvoxels*OPTIM.nconditions;
indiv_values=assigned_values(1+ig+iv:end);

%Assign global values to P
for i=1:numel(OPTIM.global_vars)
    P.(OPTIM.global_vars{i}) = global_values(i);
end

for voxel=1:OPTIM.nvoxels
    
        for j=1:numel(OPTIM.vary_voxel)
            v=OPTIM.vary_voxel{j};
            P.(v)=OPTIM.varyval_voxel(j,voxel);
        end
        
        for i=1:numel(OPTIM.voxel_vars)
            P.(OPTIM.voxel_vars{i}) = voxel_values(sub2ind([OPTIM.nvoxels, numel(OPTIM.voxel_vars)], voxel,i));
        end


        for icond=1:OPTIM.nconditions

            for jj=1:numel(OPTIM.vary_indiv)
            v=OPTIM.vary_indiv{jj};
            P.(v)=OPTIM.varyval_indiv(jj,icond);
            end
            
            for ii=1:numel(OPTIM.indiv_vars)
                P.(OPTIM.indiv_vars{ii}) = indiv_values(sub2ind([OPTIM.nvoxels, OPTIM.nconditions, numel(OPTIM.indiv_vars)], voxel,icond,ii));
            end
    
            Msim   = sim(P);

            if isfield(P,'M0corr'); Msim.zspec=Msim.zspec/P.M0corr; end
    
            if isfield(P,'normalized')
                if ~isempty(P.normalized)
                    Pnorm=P;
                    Pnorm.xZspec=P.normalized;
                    M0     = sim(Pnorm);
                    Msim.zspec=Msim.zspec./M0.zspec;
                end
            end
    
            %Adapt output depending on fit type
            if OPTIM.fit_type==0
                %zspec   = [zspec  , Msim.zspec];
                zspec(voxel,icond,:)=Msim.zspec;
                %xZspec  = [xZspec ,  P.xZspec];
            elseif OPTIM.fit_type==1
                [MTR,xzspecMTR]=calc_MTRasym(Msim.zspec,P.xZspec);
                %zspec   = [zspec  , MTR];
                zspec(voxel,icond,:)=MTR;
                %xZspec  = [xZspec ,  xzspecMTR];
            elseif OPTIM.fit_type==2
                [MTR,~]=calc_MTRasym(Msim.zspec,P.xZspec);
                wholedata   = [Msim.zspec;  MTR*OPTIM.rescalefactMTR];
                zspec(voxel,icond,:)=wholedata;
            end
        end
end





