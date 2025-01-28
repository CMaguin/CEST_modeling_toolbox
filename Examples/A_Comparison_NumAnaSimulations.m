%% Setting up simulation scheme & generic parameters
%close all

clear all
for n=1:2
P=init_Sim_struct();

%P=add_pool_model(P,'glucose_Zaiss');
% P=add_pool_model(P,'glucose_Zaiss');
if n>1
    P=add_pool_model(P,'creatine');
end
% P=add_pool_model(P,'Myo_inositol');
% P=add_pool_model(P,'Phosphocreatine');
% P=add_pool_model(P,'GABA');
% P=add_pool_model(P,'Taurine');
% P=add_pool_model(P,'Glutamine');
P=add_pool_model(P,'glutamate_singlepool_pHdependent');
% P=add_pool_model(P,'mAmides');
% P=add_pool_model(P,'NOE_-1.75ppm');
% P=add_pool_model(P,'NOE_-2.25ppm');
% P=add_pool_model(P,'NOE_-2.75ppm');
% P=add_pool_model(P,'NOE_-3.25ppm');
% P=add_pool_model(P,'NOE_-3.75ppm');
% P=add_pool_model(P,'MT'); P.MT_cutoff=3000;
% P.Glc=50;
% P.kex(strcmp(P.pool_names, 'Glc_1'))= 323; 
% P.kex(strcmp(P.pool_names, 'Glc_2')) = 492;
% P.kex(strcmp(P.pool_names, 'Glc_3'))= 374;
% P.kex(strcmp(P.pool_names, 'Glc_4'))= 675;
% P.T1_water=0.5;
% P.T2_water=0.3;

% P.Offsets_order='successive';

%Sequence parameters
P.xZspec=[-5:0.2:5];

P.FREQ=11.7*gamma_; %B0=11.7T
P.Zi=1.0;
P.B1=5;		 % the saturation B1 in µT

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
P.normalized=+20;

P.Glu=30;
P.Cr=30;

T1A=2.1;
T2A=0.03;

P.T1_water=T1A;
P.T2_water=T2A;


% Vary one parameter
vary='B1';
%vary='Glu';
%vary_values=1/(2*55.6*1000) *2*[32,16,8,4,2];
%T1A=linspace(1,4,15);
%vary_values=1./T1A;
%vary_values=[10,100,500,1000,2000,3000,4000];
%vary_values=[10,10.5];
vary_values=[1,3,5,7];
%Noise_s=0.001; %std is 0.005
Noise_s=0.005;

% Making simulations

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
    if n==1
        T(1)=plot(P.xZspec, sol_num.zspec,	'Color',"#A2142F",'LineWidth',1.0); hold on;
        T(2)=plot(P.xZspec, sol_ana.zspec,'s',	'Color',"#A2142F",'LineWidth',1.0); hold on;
        figure(2)
        plot(P.xZspec, Simulation.MTRasym(i,:),'Color',	"#A2142F",'LineWidth',1.0); hold on;
        plot(P.xZspec, Simulation_ana.MTRasym(i,:),'s',	'Color',"#A2142F",'LineWidth',1.0); hold on;
        xlim([0 max(P.xZspec)])
    else
        T(3)=plot(P.xZspec, sol_num.zspec,	'Color',"#77AC30",'LineWidth',1.0); hold on;
        T(4)=plot(P.xZspec, sol_ana.zspec,'s',	'Color',"#77AC30",'LineWidth',1.0); hold on;
        figure(2)
        plot(P.xZspec, Simulation.MTRasym(i,:),	'Color',"#77AC30",'LineWidth',1.0); hold on;
        plot(P.xZspec, Simulation_ana.MTRasym(i,:),'s',	'Color',"#77AC30",'LineWidth',1.0); hold on;
        xlim([0 max(P.xZspec)])
    end

    %figure(2)
    %plot(P.xZspec, Simulation_noisy.Zspec(i,:)); hold on;    
end


end

figure(1)
% l=legend(num2str(vary_values', '%.1f'),'Location','southeast');
% l.Title.String=vary;

set(gca,'Box','off')
set(gca,'Color','white')
set(gca,'FontSize',12)

legend(T,{'Solution numérique (2 compartiments)','Solution analytique (2 compartiments)','Solution numérique (3 compartiments)','Solution analytique (3 compartiments)'},'FontSize',12)
xlabel '\delta (ppm)'
ylabel 'Z'
%title 'Simulations des équations de Bloch-McConnell'
set(gca,'XDir','reverse');
figure(2)
% l=legend(num2str(vary_values', '%.1f'),'Location','southeast');
% l.Title.String='pH';
xlabel '\delta (ppm)'
ylabel 'MTR_{asym}'
%title 'W+Glu '
set(gca,'XDir','reverse');
set(gca,'Box','off')
set(gca,'Color','white')
set(gca,'FontSize',12)


%% Computing time
Nmax=15;
for n=1:Nmax
P=init_Sim_struct();

%P=add_pool_model(P,'glucose_Zaiss');
% P=add_pool_model(P,'glucose_Zaiss');
for iii=1:n
%     P=add_pool_model(P,'creatine');
    P=add_pool_model(P,'glutamate_physio');
end
% P=add_pool_model(P,'Myo_inositol');
% P=add_pool_model(P,'Phosphocreatine');
% P=add_pool_model(P,'GABA');
% P=add_pool_model(P,'Taurine');
% P=add_pool_model(P,'Glutamine');
%P=add_pool_model(P,'glutamate_singlepool_pHdependent');
% P=add_pool_model(P,'mAmides');
% P=add_pool_model(P,'NOE_-1.75ppm');
% P=add_pool_model(P,'NOE_-2.25ppm');
% P=add_pool_model(P,'NOE_-2.75ppm');
% P=add_pool_model(P,'NOE_-3.25ppm');
% P=add_pool_model(P,'NOE_-3.75ppm');
% P=add_pool_model(P,'MT'); P.MT_cutoff=3000;
% P.Glc=50;
% P.kex(strcmp(P.pool_names, 'Glc_1'))= 323; 
% P.kex(strcmp(P.pool_names, 'Glc_2')) = 492;
% P.kex(strcmp(P.pool_names, 'Glc_3'))= 374;
% P.kex(strcmp(P.pool_names, 'Glc_4'))= 675;
% P.T1_water=0.5;
% P.T2_water=0.3;

% P.Offsets_order='successive';

%Sequence parameters
P.xZspec=[-5:0.2:5];

P.FREQ=11.7*gamma_; %B0=11.7T
P.Zi=1.0;
P.B1=5;		 % the saturation B1 in µT

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
P.normalized=+20;

P.Glu=30;
P.Cr=30;

T1A=2.1;
T2A=0.03;

P.T1_water=T1A;
P.T2_water=T2A;


% Vary one parameter
vary='B1';
%vary='Glu';
%vary_values=1/(2*55.6*1000) *2*[32,16,8,4,2];
%T1A=linspace(1,4,15);
%vary_values=1./T1A;
%vary_values=[10,100,500,1000,2000,3000,4000];
%vary_values=[10,10.5];
vary_values=[3];
%Noise_s=0.001; %std is 0.005
Noise_s=0.005;

% Making simulations

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
    
    tic
    sol_ana =analytic_simulation(P);
    Computing_time_ana(n)=toc;
    
    tic
    sol_num =numeric_simulation(P);
    Computing_time_num(n)=toc;

    Difference_numana(n)=mean(abs(sol_ana.zspec-sol_num.zspec));% sqrt(mean((sol_ana.zspec-sol_num.zspec).^2));
    
    Simulation.Zspec(i, :)=sol_num.zspec;
    Simulation.MTRasym(i,:)=sol_num.zspec(end:-1:1)-sol_num.zspec;
    Simulation_ana.Zspec(i, :)=sol_ana.zspec;
    Simulation_ana.MTRasym(i,:)=sol_ana.zspec(end:-1:1)-sol_ana.zspec;
%      %Plotting
%     figure(1)
%     if n==1
%         T(1)=plot(P.xZspec, sol_num.zspec,'r'); hold on;
%         T(2)=plot(P.xZspec, sol_ana.zspec,'b+'); hold on;
%         figure(2)
%         plot(P.xZspec, Simulation.MTRasym(i,:),'r'); hold on;
%         plot(P.xZspec, Simulation_ana.MTRasym(i,:),'b+'); hold on;
%         xlim([0 max(P.xZspec)])
%     else
%         T(3)=plot(P.xZspec, sol_num.zspec,'g'); hold on;
%         T(4)=plot(P.xZspec, sol_ana.zspec,'bs'); hold on;
%         figure(2)
%         plot(P.xZspec, Simulation.MTRasym(i,:),'g'); hold on;
%         plot(P.xZspec, Simulation_ana.MTRasym(i,:),'bs'); hold on;
%         xlim([0 max(P.xZspec)])
%     end

    %figure(2)
    %plot(P.xZspec, Simulation_noisy.Zspec(i,:)); hold on;    
end


end

%
figure()

%yyaxis left
hold on
plot(1:Nmax,Computing_time_num,'r','LineWidth',1.5)
plot(1:Nmax,Computing_time_ana,'b','LineWidth',1.5)
set(gca,'YScale','log')

ylabel 'Temps de calcul (s)'
xlabel 'Nombre de compartiments dans le modèle'

ylim([9*10^(-4), 1])
set(gca,'Box','off')
set(gca,'Color','white')
set(gca,'FontSize',12)
legend('Solution numérique','Solution analytique')

% yyaxis right
% plot(1:Nmax,Difference_numana,'k')

