%% Setting up simulation scheme & generic parameters
%close all

% clear all
P=init_Sim_struct();

% 
% P=add_pool_model(P,'glucose_Zaiss');
% P=add_pool_model(P,'creatine');
% % P=add_pool_model(P,'Myo_inositol');
% P=add_pool_model(P,'Phosphocreatine');
% P=add_pool_model(P,'GABA');
% P=add_pool_model(P,'Taurine');
% P=add_pool_model(P,'Glutamine');
% P=add_pool_model(P,'glutamate_singlepool_pHdependent');
% P=add_pool_model(P,'mAmides');
% % P=add_pool_model(P,'NOE_-1.75ppm');
% % P=add_pool_model(P,'NOE_-2.25ppm');
% % P=add_pool_model(P,'NOE_-2.75ppm');
% % P=add_pool_model(P,'NOE_-3.25ppm');
% P=add_pool_model(P,'NOE_-3.75ppm');
P=add_pool_model(P,'MT');

P.T1_water=2;
P.T2_water=0.03;


P.MT_lineshape='SuperLorentzian';
P.MT_cutoff=0;

P.H_T=0.05*P.H_water;
P.T2_MT=10*10^(-6);
P.dw_MT=-2.34;

% P.Offsets_order='successive';

%Sequence parameters
%P.xZspec=[-100;-50;-30;-20;-10;-9;-8;-7;-6.50000000000000;-6;-5.50000000000000;-5;-4.50000000000000;-4;-3.50000000000000;-3;-2.50000000000000;-2;-1.50000000000000;-1;-0.500000000000000;0;0.500000000000000;1;1.50000000000000;2;2.50000000000000;3;3.50000000000000;4;4.50000000000000;5;5.50000000000000;6;6.50000000000000;7;8;9;10;20;30;50;100];
P.xZspec=-20:0.1:20;

P.FREQ=11.7*gamma_; %B0=11.7T
P.Zi=1.0;
P.B1=5;		 % the saturation B1 in µT

P.pulsed=1;
if P.pulsed
    P.shape = 'block';                  % choose 'block' for cw saturation (cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4)
    P.n     = 10* 1;                       % number of saturation pulses
    P.tp    = 0.1;                      % saturation time per pulse in s
    P.td=0.00001; %10µs
    P.td=0;
else
    P.tp=1;
    P.n=1;
    P.td=0;
    P.shape='block';
end
P.TR=5; %s
P.Trec=P.TR-(P.tp+P.td)*P.n; %s
%P.normalized=-100;


%% Vary one parameter
vary='MT_cutoff';
vary_values=P.FREQ*2*pi*[1,1,1,2,4]; %1,2,5ppm cutoff

%Noise_s=0.001; %std is 0.005
Noise_s=0.0005;

% Making simulations

%Set up a structure to save simulated Z spectrum
Simulation.Sim=P;
Simulation.vary=vary;
Simulation.vary_values=vary_values;
Simulation.Zspec=zeros(numel(vary_values), numel(P.xZspec));
Simulation.MTRasym=zeros(numel(vary_values), numel(P.xZspec));

%Add noise and set up a structure to save simulated noisy Z spectrum
Simulation_noisy.Sim=P;
Simulation_noisy.vary=vary;
Simulation_noisy.vary_values=vary_values;
Simulation_noisy.Zspec=zeros(numel(vary_values), numel(P.xZspec));
Simulation_noisy.MTRasym=zeros(numel(vary_values), numel(P.xZspec));
Simulations_noisy.s_noise=Noise_s;

Colors={"#0072BD",	"#D95319",	"#EDB120","#7E2F8E","#77AC30"};

for i=1:5
    if i>2
        P.(vary)=vary_values(i);
        P.MT_lineshape='SuperLorentzian';
    elseif i==2
        P.MT_lineshape='Lorentzian';
    else
        P.MT_lineshape='Gaussian';
    end
    P=update_P(P);
    if P.analytic==1
        sol =analytic_simulation(P);
    else
        sol =numerci_simulation(P);
    end
    Simulation.Zspec(i, :)=sol.zspec;
    Simulation.MTRasym(i,:)=sol.zspec(end:-1:1)-sol.zspec;
    Simulation_noisy.Zspec(i,:)=ricernd(Simulation.Zspec(i,:), Noise_s);
    Simulation_noisy.MTRasym(i,:)=Simulation_noisy.Zspec(i,end:-1:1)-Simulation_noisy.Zspec(i,:);
    
    %Plotting
    figure(1)
    T(i)=plot(P.xZspec, sol.zspec,'LineWidth',1.5,'Color',Colors{i}); hold on;
    set(gca,'XDir','reverse');
    figure(2)
    yyaxis left
    plot(P.xZspec, sol.zspec,'LineWidth',1.5,'Color',Colors{i},'LineStyle','-','Marker','none'); hold on;
    posppm=(P.xZspec<=0);
    MTR=sol.zspec(end:-1:1)-sol.zspec;
    xlabel '\delta (ppm)'
    ylabel 'Z'
    yyaxis right
    xlabel '\delta (ppm)'
    ylabel 'MTR_{asym} (%)'
    plot(P.xZspec(posppm), 100*MTR(posppm),'LineWidth',1.5,'Color',Colors{i},'LineStyle','-','Marker','none'); hold on;
    set(gca,'XDir','reverse');
    ylim([0 10])
    xlim([-5 5])
    figure(3)
    plot(P.xZspec, sol.S.Rex_MT,'LineWidth',1.5,'SeriesIndex',i); hold on;
    set(gca,'XDir','reverse');
    xlabel '\delta (ppm)'
    ylabel 'R_{ex}^{MT} (Hz)'
end
figure(1)
plot([-5,5,5,-5,-5],[0.01,0.01,1,1,0.01],'k--') % Zoom box
l=legend(T,{'Gaussienne','Lorentzienne','Super-Lorentzienne \newline (coupure \pm 1 ppm)','Super-Lorentzienne \newline (coupure \pm 2 ppm)','Super-Lorentzienne \newline (coupure \pm 5 ppm)'},'Location','southeast');
l.Title.String='Modélisations MT';
set(gca,'Box','off')
ylabel 'Z'
xlabel '\delta (ppm)'
set(gca,'XDir','reverse');
figure(2)
set(gca,'Box','off')
xlabel '\delta (ppm)'
set(gca,'XDir','reverse');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

