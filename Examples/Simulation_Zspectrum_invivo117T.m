%% Setting up simulation scheme & generic parameters
close all

clear all
P=init_Sim_struct();

% 
P=add_pool_model(P,'creatine');
P=add_pool_model(P,'Myo_inositol');
P=add_pool_model(P,'Phosphocreatine');
P=add_pool_model(P,'GABA');
P=add_pool_model(P,'Taurine');
P=add_pool_model(P,'Glutamine');
P=add_pool_model(P,'glutamate_singlepool_pHdependent');
P=add_pool_model(P,'mAmides');
P=add_pool_model(P,'NOE_-3.5ppm');
% P=add_pool_model(P,'NOE_-1.6ppm');
P=add_pool_model(P,'MT');


% P.Offsets_order='successive';

%Sequence parameters
%P.xZspec=[-100;-50;-30;-20;-10;-9;-8;-7;-6.50000000000000;-6;-5.50000000000000;-5;-4.50000000000000;-4;-3.50000000000000;-3;-2.50000000000000;-2;-1.50000000000000;-1;-0.500000000000000;0;0.500000000000000;1;1.50000000000000;2;2.50000000000000;3;3.50000000000000;4;4.50000000000000;5;5.50000000000000;6;6.50000000000000;7;8;9;10;20;30;50;100];
P.xZspec=-5:0.1:5;

P.FREQ=11.7*gamma_; %B0=11.7T
P.Zi=1.0;
P.B1=1;		 % the saturation B1 in µT

P.pulsed=0;
if P.pulsed
    P.shape = 'block';                  % choose 'block' for cw saturation (cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4)
    P.n     = 10* 1;                       % number of saturation pulses
    P.tp    = 0.1;                      % saturation time per pulse in s
    P.td=0.00001; %10µs
else
    P.tp=1;
    P.n=1;
    P.td=0;
    P.shape='block';
end
P.TR=5; %s
P.Trec=P.TR-(P.tp+P.td)*P.n; %s
P.normalized=-100;


%% Vary one parameter
vary='B1';
vary_values=[1,3,5,7];

%Noise_s=0.001; 
Noise_s=0.001;

%% Making simulations

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

for i=1:numel(vary_values)
    P.(vary)=vary_values(i);
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
    
    c=diverging_colors_2(round(vary_values(i)),10);

    figure(1)
    yyaxis left
    T(i)=plot(P.xZspec, Simulation.Zspec(i, :),'LineWidth',2.0,'LineStyle','-','Marker','none','Color',c); hold on;
    yyaxis right
    posppm=(P.xZspec>=0);
    plot(P.xZspec(posppm), 100*Simulation.MTRasym(i,posppm),'LineWidth',1.5,'LineStyle','--','Marker','none','Color',c); hold on;
    label_list{i}=strcat(num2str(vary_values(i),'%.0f'), 'µT');

    figure(2)
    yyaxis left
    T2(i)=plot(P.xZspec, Simulation_noisy.Zspec(i,:),'LineWidth',2.0,'LineStyle','-','Marker','none','Color',c); hold on;
    yyaxis right
    posppm=(P.xZspec>=0);
    plot(P.xZspec(posppm), 100*Simulation_noisy.MTRasym(i,posppm),'LineWidth',1.5,'LineStyle','--','Marker','none','Color',c); hold on;
end
figure(1)
l=legend(T,label_list,'Location','southeast');
l.Title.String='B_1 power';
xlabel 'Offset (ppm)'
yyaxis left
ylim([0 1])
yticks(0:0.1:1)
set(gca,'YColor','k')
ylabel 'Z'
set(gca,'XDir','reverse');
yyaxis right
ylim([-20 40])
yticks(0:5:15)
set(gca,'YColor','k')
title 'Z-spectra in vivo'
figure(2)
l2=legend(T2,label_list,'Location','southeast');
l2.Title.String='B_1 power';
xlabel 'Offset (ppm)'
yyaxis left
ylim([0 1])
yticks(0:0.1:1)
set(gca,'YColor','k')
ylabel 'Z'
set(gca,'XDir','reverse');
yyaxis right
ylim([-20 40])
yticks(0:5:15)
set(gca,'YColor','k')
title 'Z-spectra in vivo but with some noise'
set(gca,'XDir','reverse');


   