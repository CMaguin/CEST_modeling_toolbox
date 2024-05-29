function zspec=compute_ODE_lipoGlu(P, OPTIM, assigned_values,xdata)

if nargin>3
    P.xZspec=xdata;
end

% if isfield(P, 'sim_type')
%     if P.sim_type==0 
%         sim=@(P) numeric_simulation(P); 
%     else
%         sim=@(P) analytic_simulation(P); 
%     end
% else
       % sim=@(P) analytic_simulation(P); 
% end
    

sim=@(P) ODE_simulation_lipoGlu(P);

if iscolumn(P.xZspec)
    P.xZspec=P.xZspec';
end

multiple=0;
if isfield(OPTIM,'vary')
   if ~strcmp(OPTIM.vary,'None')
       multiple=1;
       val=OPTIM.vary_val;
   end
end


%Assign values to P
for i=1:numel(OPTIM.dep_vars)
    P.(OPTIM.dep_vars{i}) = assigned_values(i);
end


if multiple
    zspec=[];
    xZspec=[];
    for ii=1:numel(val(1,:))
        for jj = 1:numel(OPTIM.vary) % this realizes multiple varying parameters
            variable=OPTIM.vary{jj};
            P.(variable) = val(jj,ii);
        end
        
        
            Msim=sim(P);

            if isfield(P,'M0corr'); Msim.zspec=Msim.zspec/P.M0corr; end

            if isfield(P,'normalized')
                if ~isempty(P.normalized)
                    Pnorm=P;
                    Pnorm.xZspec=P.normalized;
                    M0     = ODE_simulation_lipoGlu(Pnorm);
                    Msim.zspec=Msim.zspec./M0.zspec;
                end
            end

            %Adapt output depending on fit type
            if OPTIM.fit_type==0
                zspec   = [zspec  , Msim.zspec];
                xZspec  = [xZspec ,  P.xZspec];
            elseif OPTIM.fit_type==1
                [MTR,xzspecMTR]=calc_MTRasym(Msim.zspec,P.xZspec);
                zspec   = [zspec  , MTR];
                xZspec  = [xZspec ,  xzspecMTR];
            elseif OPTIM.fit_type==2
                [MTR,~]=calc_MTRasym(Msim.zspec,P.xZspec);
                wholedata   = [Msim.zspec;  MTR*OPTIM.rescalefactMTR];
                zspec   = [zspec  , wholedata];
            end
    end
    zspec=zspec';
else
    Msim=sim(P);
    if isfield(P,'M0corr'); Msim.zspec=Msim.zspec/P.M0corr; end
    if isfield(P,'normalized')
        if ~isempty(P.normalized)
            Pnorm=P;
            Pnorm.xZspec=P.normalized;
            M0     = ODE_simulation_lipoGlu(Pnorm);
            Msim.zspec=Msim.zspec./M0.zspec;
        end
    end
    %Adapt output depending on fit type
    if OPTIM.fit_type==0
        zspec   = Msim.zspec;
        xZspec  = P.xZspec;
    elseif OPTIM.fit_type==1
        [MTR,xzspecMTR]=calc_MTRasym(Msim.zspec,P.xZspec);
        zspec   = MTR;
        xZspec  = xzspecMTR;
    elseif OPTIM.fit_type==2
        [MTR,~]=calc_MTRasym(Msim.zspec,P.xZspec);
        wholedata   = [Msim.zspec;  MTR*OPTIM.rescalefactMTR];
        zspec   = wholedata;
    end

end;


