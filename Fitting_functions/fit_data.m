function [Results,Fit, Pnew]=fit_data(P,Z,OPTIM,options)

if nargin<4 %default fit options
    options = optimset('TolFun',1e-18,'TolX', 1e-18,'MaxIter',max(400, 200*numel(OPTIM.start)),'MaxFunEvals',5000*numel(OPTIM.start),'Display','on');
end


if OPTIM.fit_type==0
    ydata=Z;
    xdata=P.xZspec;
elseif OPTIM.fit_type==1
    [ydata,xdata]=calc_MTRasym(Z,P.xZspec);
elseif OPTIM.fit_type==2
    [MTR,xMTR]=calc_MTRasym(Z,P.xZspec);
    if numel(size(Z))>2
    ydata=cat(3,Z, MTR*OPTIM.rescalefactMTR);
    else
    ydata=[Z MTR*OPTIM.rescalefactMTR];
    end
    xdata=[P.xZspec xMTR];
end

ydata=reshape(ydata, OPTIM.nvoxels,OPTIM.nconditions, []);

if isfield(OPTIM, 'multivar')
if OPTIM.multivar
    f = @(x,xdata) compute_analytical_multivar(P, OPTIM, x);
end
else
    f = @(x,xdata) compute_analytical(P, OPTIM, x);
end


[x, resnorm, RES,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(f,OPTIM.start,P.xZspec,ydata,OPTIM.lowerbounds,OPTIM.upperbounds,options);

SStot = sum((ydata(:)-mean(ydata(:))).^2); 
R_sq= 1-(resnorm/SStot);
Results.fitted_param_vector=x;
Results.R_sq=R_sq;
Results.SStot=SStot;
Results.J=JACOBIAN;
Results.RSS=resnorm;

% ----Information criterions - help estimate which model is the best----
% (assess how many degrees of freedom needed while "avoiding" overfitting)
%Akaike information criterion
sigma=resnorm/numel(ydata);
AIC=2*numel(x)+numel(ydata)*log(sigma) + numel(ydata)*(1+log(2*pi)); %last term is not important but let's add it for the sake of math
Results.AIC=AIC;
% Corrected Akaike for sample samples
try
    cAIC=AIC + 2*numel(x)*(numel(x)+1)/(numel(ydata)-numel(x)-1);
catch
    cAIC=NaN;
end
Results.cAIC=cAIC;
%Bayesian Information Criterion
BIC=numel(ydata)*log(sigma)+log(numel(ydata))*numel(x) + numel(ydata)*(1+log(2*pi));% last term is not important but let's be consistant
Results.BIC=BIC;

%Calculate error
[ci, varb, corrb, varinf] = nlparci_custom(x,resnorm,JACOBIAN,0.5);
ci(:,1) = ci(:,1)-x;
ci(:,2) = ci(:,2)-x;

%Assign results to Pnew
Pnew=P;
for i=1:numel(OPTIM.dep_vars)
    Pnew.(OPTIM.dep_vars{i})=x(i);
end
Pnew=update_P(Pnew);

%Results structure
for i=1:numel(OPTIM.dep_vars)
    Results.(OPTIM.dep_vars{i})=x(i);
end
Results.ci=ci;
Results.resnorm=resnorm;
if isfield(OPTIM, 'multivar')
if OPTIM.multivar
    Results.voxel_vars=OPTIM.voxel_vars;
    Results.voxel_val=x(1+numel(OPTIM.global_vars):numel(OPTIM.global_vars)+numel(OPTIM.voxel_vars)*OPTIM.nvoxels);
    Results.ci_voxel=ci(1+numel(OPTIM.global_vars):numel(OPTIM.global_vars)+numel(OPTIM.voxel_vars)*OPTIM.nvoxels,:);
    Results.voxel_val=reshape(Results.voxel_val, OPTIM.nvoxels,numel(OPTIM.voxel_vars));
    Results.ci_voxel=reshape(Results.ci_voxel,OPTIM.nvoxels,numel(OPTIM.voxel_vars),2);
    Results.indiv_vars=OPTIM.indiv_vars;
    Results.indiv_val=x(1+numel(OPTIM.global_vars)+numel(OPTIM.voxel_vars)*OPTIM.nvoxels:end);
    Results.ci_indiv=ci(1+numel(OPTIM.global_vars)+numel(OPTIM.voxel_vars)*OPTIM.nvoxels:end,:);
    Results.ci_indiv=reshape(Results.ci_indiv,OPTIM.nvoxels,OPTIM.nconditions,numel(OPTIM.indiv_vars),2);
    Results.indiv_val=reshape(Results.indiv_val,OPTIM.nvoxels,OPTIM.nconditions,numel(OPTIM.indiv_vars));
end
end



%Create fit structure
Fit.x=P.xZspec;
Fit.fit=f(x,P.xZspec);
Fit.Rsq=R_sq;
try
    densexxx=min(P.xZspec):0.1:max(P.xZspec);
    noffofsets_fit=numel(densexxx);
    Fit.densex=densexxx;

    
    if isfield(OPTIM, 'multivar')
    if OPTIM.multivar
        f = @(x,xdata) compute_analytical_multivar(P, OPTIM, x,densexxx);
    end
    else
        f = @(x,xdata) compute_analytical(P, OPTIM, x,densexxx);
    end
    Fit.fitdense=f(x,densexxx);
catch
    disp('Unable to do return dense fit')
end
%Fit.fit=f(x,densexxx);
% Fit.fit=interp1(P.xZspec, Fit.fit, densexxx);

% if isfield(OPTIM, 'multivar')
% if OPTIM.multivar
%     f = @(x,xdata) compute_analytical_multivar(P, OPTIM, x);
% end
% else
%     f = @(x,xdata) compute_analytical(P, OPTIM, x);
% end
% if  OPTIM.fit_type~=2
%     Fit.fit=f(x,densexxx);
% else 
%     whole_fit=f(x,densexxx);
%     whole_fit=reshape(whole_fit,OPTIM.nvoxels, OPTIM.nconditions, []);
%     Fit.fitZ=whole_fit(:,:,1:floor(noffofsets_fit));
%     Fit.fitMTR=whole_fit(:,:,floor(noffofsets_fit)+1:end)/OPTIM.rescalefactMTR;
% end
Fit.output=OUTPUT;
