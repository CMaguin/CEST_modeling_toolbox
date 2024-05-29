%% Setting up simulation scheme & generic parameters
%close all

clear all
P=init_Sim_struct();

P=add_pool_model(P,'glutamate_physio');
%P=add_pool_model(P,'creatine');
% P=add_pool_model(P,'Glutamine');
% P=add_pool_model(P,'Taurine');
% P=add_pool_model(P,'GABA');
%P=add_pool_model(P,'MT');

P.MT=0;
%P.MT_lineshape='SuperLorentzian';

T1A=2.2;
T2A=0.05;

P.T1_water=T1A;
P.T2_water=T2A;

% P.Offsets_order='successive';

%Sequence parameters
P.xZspec=[-5:0.1:5];

P.FREQ=11.7*gamma_; %B0=11.7T
P.Zi=1.0;
P.B1=1;		 % the saturation B1 in µT

P.pulsed=0;
if P.pulsed
    P.shape = 'block';                  % choose 'block' for cw saturation (cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4)
    P.n     = 10;                       % number of saturation pulses
    P.tp    = 0.1;                      % saturation time per pulse in s
    P.td=0.01;
else
    P.tp=1;
    P.n=1;
    P.td=0;
    P.shape='block';
end
P.TR=5; %s
P.Trec=P.TR-(P.tp+P.td)*P.n; %s
P.normalized=-20;

%P.sim_type=0;

P.kex_Glu=7400;
P.Glu=30;

P.CALC=calc_alt_var(P);

%% Vary one parameter
vary='B1';
%vary_values=1/(2*55.6*1000) *2*[32,16,8,4,2];
%T1A=linspace(1,4,15);
%vary_values=1./T1A;
vary_values=1:6;



%% Making simulations

%Set up a structure to save simulated Z spectrum
Simulation.Sim=P;
Simulation.vary=vary;
Simulation.vary_values=vary_values;
Simulation.Zspec=zeros(numel(vary_values), numel(P.xZspec));
Simulation.MTRasym=zeros(numel(vary_values), numel(P.xZspec));
Simulation_ana.Zspec=zeros(numel(vary_values), numel(P.xZspec));
Simulation_ana.MTRasym=zeros(numel(vary_values), numel(P.xZspec));


for i=1:numel(vary_values)
    P.(vary)=vary_values(i);
    P=update_P(P);
    
    sol_ana =analytic_simulation(P);
    sol_num =numeric_simulation(P);
    
    Simulation.Zspec(i, :)=sol_num.zspec;
    Simulation.MTRasym(i,:)=sol_num.zspec(end:-1:1)-sol_num.zspec;
    Simulation_ana.Zspec(i, :)=sol_ana.zspec;
    Simulation_ana.MTRasym(i,:)=sol_ana.zspec(end:-1:1)-sol_ana.zspec;
     %Plotting
    figure(1)
    plot(P.xZspec, sol_num.zspec,'r'); hold on;
    plot(P.xZspec, sol_ana.zspec,'b+'); hold on;
    figure(2)
    plot(P.xZspec, Simulation.MTRasym(i,:),'r'); hold on;
    plot(P.xZspec, Simulation_ana.MTRasym(i,:),'b+'); hold on;
    xlim([0 max(P.xZspec)])
    %figure(2)
    %plot(P.xZspec, Simulation_noisy.Zspec(i,:)); hold on;    
end
figure(1)
% l=legend(num2str(vary_values', '%.1f'),'Location','southeast');
% l.Title.String=vary;

ax=gca;
ax.YColor='k';
ax.color='white';

legend('Numeric solution','Analytic solution')
xlabel '\delta (ppm)'
ylabel 'Z'
title 'W+Glu'
set(gca,'XDir','reverse');
figure(2)
% l=legend(num2str(vary_values', '%.1f'),'Location','southeast');
% l.Title.String='pH';
xlabel 'Offset (ppm)'
ylabel 'MTR_{asym}'
title 'W+Glu '
set(gca,'XDir','reverse');

%% Setting up simulation scheme & generic parameters
%close all
P=init_Sim();
%Sequence parameters
P.analytic=0;
P.Rex_sol='Lorentz';
P.xZspec=[-10:0.2:10];

P.FREQ=11.7*gamma_; %B0=11.7T
P.Zi=1.0;
P.B1=1;		 % the saturation B1 in µT
P.spoilf=1;

P.dummies=1;

P.pulsed=1;
if P.pulsed
    P.shape = 'block';                  % choose 'block' for cw saturation (cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4)
    P.n     = 10;                       % number of saturation pulses
    P.tp    = 0.1;                      % saturation time per pulse in s
    P.td=0.01;
    P.DC=0.9091;
else
    P.tp=1;
    P.n=1;
    P.td=0;
    P.shape='block';
end
P.TR=5; %s
P.Trec=P.TR-(P.tp+P.td)*P.n; %s
P.normalized=-20;

P.dummies=0;

%Set up nbr of pools (excluding water & MT)
P.n_cest_pool=2;

%Choose if MT and parameters
P.MT=1;
P.MT_sol_type='Rex_MT';
P.MT_shape='SuperLorentzian';
P.MT_lineshape='SuperLorentzian';

%Set up your different pools & physiological parameters
P.CESTparameters='custom'; %generic if you want values as described in fct getSim, 'custom' if you want to set up parameters as you want
P.tissue='GM';
P.CESTagent='Creatine';
P.n_cest_pools=1;

if strcmp(P.CESTparameters, 'generic')
    P=getSim(P); 
elseif strcmp(P.CESTparameters, 'custom') % change here if you want custom parameters
    P=init_Sim_param(P);
    P.n_cest_pool   = 2;
    %Water parameters
    P.dwA=0;
    P.R1A=1/T1A;
    P.R2A=1/T2A;

    %Pool B parameters
    P.fB=3*10/(2*55.6*1000); % fb=nbr of protons of CEST * C(mM)/(nbr of protons of water * Cwater in 1L *1000mL)
    P.kBA=7738;
    P.R2B=66.6667;
    P.dwB=3.0;
    P.kAB=P.fB*P.kBA;
%     P.R2B=1/T2A;
    P.R1B=1/T1A;
    
    %Pool D parameters
    P.fD=4*5/(2*55.6*1000); % fb=nbr of protons of CEST * C(mM)/(nbr of protons of water * Cwater in 1L *1000mL)
    P.kDA=810;
    P.R2D=128.2051;
    P.dwD=2.0;
    P.kAD=P.fD*P.kDA;
%     P.R2D=1/T2A;
    P.R1D=1/T1A;

    %MT parameters
    P.MT=1;
    P.fC=0.0495;
    P.dwC=-2.34;
    P.R2C=1e5;
    P.R1C=1;
    P.kCA=40;
    P.kAC=P.fC*P.kCA;

    P=calc_Hconc(P);
    P=calc_fH(P);
    P=calc_backexch_rates(P);
end

P=calc_Hconc(P);
P=calc_backexch_rates(P);



%% Making simulations

%Set up a structure to save simulated Z spectrum
Simulationold.Sim=P;
Simulationold.vary=vary;
Simulationold.vary_values=vary_values;
Simulationold.Zspec=zeros(numel(vary_values), numel(P.xZspec));
Simulationold.MTRasym=zeros(numel(vary_values), numel(P.xZspec));



for i=1:numel(vary_values)
    P.(vary)=vary_values(i);
    if P.analytic==1
        sol =ANALYTIC_SIM(P);
    else
        sol =NUMERIC_SIM(P);
    end
    Simulationold.Zspec(i, :)=sol.zspec;
    Simulationold.MTRasym(i,:)=sol.zspec(end:-1:1)-sol.zspec;
    %Plotting
    figure(1)
    plot(P.xZspec, sol.zspec,'b'); hold on;
    figure(2)
    plot(P.xZspec, Simulationold.MTRasym(i,:)); hold on;
    %figure(2)
    %plot(P.xZspec, Simulation_noisy.Zspec(i,:)); hold on;    
end
figure(1)
l=legend(num2str(vary_values', '%.1f'),'Location','southeast');
l.Title.String='B1 (µT)';
xlabel 'Offset (ppm)'
ylabel 'Z'
title 'W+Glu+Glc ([Glu]=10mM, [Glc]=3mM)'
figure(2)
l=legend(num2str(vary_values', '%.1f'),'Location','southeast');
l.Title.String='B1 (µT)';
xlabel 'Offset (ppm)'
ylabel 'MTR_{asym}'
title 'W+Glu+Glc ([Glu]=10mM, [Glc]=3mM)'
