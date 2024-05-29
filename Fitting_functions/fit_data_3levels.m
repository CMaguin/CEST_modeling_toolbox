function [Results,Fit, Pnew]=fit_data_3levels(P,Z,OPTIM)

global Results_indiv
global Results_voxels

OPTIM.dep_vars=OPTIM.global_vars;
OPTIM.start=OPTIM.global_start;
OPTIM.lowerbounds=OPTIM.global_lb;
OPTIM.upperbounds=OPTIM.global_ub;

s=size(Z);
if s(1)~=OPTIM.nvoxels || s(2)~=OPTIM.nconditions || s(3)~=numel(P.xZspec)
    warning('Z is not reshaped properly (need nvoxels x nconditions x noffsets');
    Z=reshape(Z,[OPTIM.nvoxels, OPTIM.nconditions, numel(P.xZspec)]);
end

if OPTIM.fit_type==0
    ydata=Z;
    xdata=P.xZspec;
elseif OPTIM.fit_type==1
    [ydata,xdata]=calc_MTRasym(Z,P.xZspec);
elseif OPTIM.fit_type==2
    [MTR,xMTR]=calc_MTRasym(Z,P.xZspec);
    ydata=cat(3,Z, MTR*OPTIM.rescalefactMTR);
    xdata=[P.xZspec xMTR];
end

f=@(x, xdata) compute_for_this_global(P,OPTIM,ydata,x) ;
options = optimset('TolFun',1e-14,'TolX', 1e-14,'MaxIter',400,'Display','on');

[x, resnorm, RES,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(f,OPTIM.start,P.xZspec,ydata(:),OPTIM.lowerbounds,OPTIM.upperbounds,options);
%Calculate error
[ci, varb, corrb, varinf] = nlparci_custom(x,resnorm,JACOBIAN,0.5);
ci(:,1) = ci(:,1)-x';
ci(:,2) = ci(:,2)-x';

%Assign results to Pnew
Pnew=P;
for i=1:numel(OPTIM.dep_vars)
    Pnew.(OPTIM.dep_vars{i})=x(i);
end

%Results structure
for i=1:numel(OPTIM.dep_vars)
    Results.(OPTIM.dep_vars{i})=x(i);
end
Results.ci=ci;
Results.resnorm=resnorm;
Results.indiv=Results_indiv;
Results.voxel=Results_voxels;

%Create fit structure
noffofsets_fit=201;
densexxx=linspace(min(P.xZspec), max(P.xZspec), noffofsets_fit)';
P.xZspec=densexxx;
Fit.x=densexxx;
if  OPTIM.fit_type~=2
    Fit.fit=f(x,densexxx);
else 
    whole_fit=f(x,densexxx);
    Fit.fitZ=whole_fit(:,:,1:floor(noffofsets_fit));
    Fit.fitMTR=whole_fit(:,:,floor(noffofsets_fit)+1:end)/OPTIM.rescalefactMTR;
end
Fit.output=OUTPUT;
end

function res=compute_for_this_global(P,OPTIM,ydata,x_global)
    global Results_voxel
    for i=1:numel(OPTIM.dep_vars)
        P.(OPTIM.dep_vars{i}) = x_global(i);
    end
    [res, Results_voxel]=all_voxels_fit(P,ydata,OPTIM);
    res=res(:);
end

function [AllvoxelsZfit, Results_voxel]=all_voxels_fit(P,ydata,OPTIM)
    OPTIM.dep_vars=OPTIM.voxel_vars;
    Results_voxel.vars=OPTIM.dep_vars;
    for voxel=1:numel(OPTIM.nvoxels)
        ydata_c=squeeze(ydata(voxel,:,:));
        %assigning right values to P for this condition
        for j=1:numel(OPTIM.vary_voxel)
            v=OPTIM.vary_voxel{j};
            P.(v)=OPTIM.varyval_voxel(j);
        end
        
        OPTIM.start=zeros(numel(OPTIM.dep_vars),1); OPTIM.lowerbounds=zeros(numel(OPTIM.dep_vars),1); OPTIM.upperbounds=zeros(numel(OPTIM.dep_vars),1);
        for k=1:numel(OPTIM.dep_vars)
            OPTIM.start(k)=eval(OPTIM.voxel_start{k});
            OPTIM.lowerbounds(k)=eval(OPTIM.voxel_lb{k});
            OPTIM.upperbounds(k)=eval(OPTIM.voxel_ub{k});
        end
    
        f = @(x,xdata) compute_for_this_voxel(P,OPTIM,ydata_c,voxel,x);
    
        options = optimset('TolFun',1e-14,'TolX', 1e-14,'MaxIter',200,'Display','off');
    
        [x, resnorm, RES,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(f,OPTIM.start,P.xZspec,ydata_c(:),OPTIM.lowerbounds,OPTIM.upperbounds,options);
        
        Results_voxel.values(voxel,:)=x;
        AllvoxelsZfit(voxel,:,:)=f(x,P.xZspec);
    end
end

function res=compute_for_this_voxel(P,OPTIM,ydata,voxelnbr,x_voxel)
    global Results_indiv
    for i=1:numel(OPTIM.dep_vars)
        P.(OPTIM.dep_vars{i}) = x_voxel(i);
    end
    [res, Results_indiv]=voxel_fit(P,OPTIM,ydata); 
    res=res(:);
end

function [allindivZ, Results_indiv]=voxel_fit(P,OPTIM,ydata)
    OPTIM.dep_vars=OPTIM.indiv_vars;
    Results_indiv.vars=OPTIM.dep_vars;
    for icond=1:OPTIM.nconditions

        ydata_c=squeeze(ydata(icond,:));

        %assigning right values to P for this condition
        for j=1:numel(OPTIM.vary_indiv)
            v=OPTIM.vary_indiv{j};
            P.(v)=OPTIM.varyval_indiv(j);
        end
        
        OPTIM.start=zeros(numel(OPTIM.dep_vars),1); OPTIM.lowerbounds=zeros(numel(OPTIM.dep_vars),1); OPTIM.upperbounds=zeros(numel(OPTIM.dep_vars),1);
        for k=1:numel(OPTIM.dep_vars)
            OPTIM.start(k)=eval(OPTIM.indiv_start{k});
            OPTIM.lowerbounds(k)=eval(OPTIM.indiv_lb{k});
            OPTIM.upperbounds(k)=eval(OPTIM.indiv_ub{k});
        end

        OPTIM.vary={};
        OPTIM.varyval=[];

        f = @(x,xdata) compute_analytical(P, OPTIM, x);
    
        options = optimset('TolFun',1e-14,'TolX', 1e-14,'MaxIter',100,'Display','off');
    
        [x, resnorm, RES,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(f,OPTIM.start,P.xZspec,ydata_c(:),OPTIM.lowerbounds,OPTIM.upperbounds,options);
        

        Results_indiv.values(icond,:)=x;
        Results_indiv.resnorm(icond)=resnorm;

        allindivZ(icond,:)=f(x,P.xZspec);
    end
end

