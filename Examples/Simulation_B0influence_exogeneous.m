%% Setting up simulation scheme & generic parameters
close all

clear all
P=init_Sim_struct();




%Sequence parameters
P.xZspec=[-100:0.05:100];

P.FREQ=11.75*gamma_; %B0=11.7T
P.Zi=1.0;
P.B1=5;		 % the saturation B1 in µT

P=add_pool_model(P,'glutamate_physio');

P.T1_water=2.0;
P.T2_water=0.03;

P.Glu=15;
P.kex=2000;
P.dw_Glu=80;
P.T2_Glu=0.01;
P.T1_Glu=2.0;

P=update_P(P);



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

P.Trec=1000;

%P.normalized=+20;

P.numeric=1;
P.analytic=0;



%% Vary one parameter
vary={'FREQ','T1_water','T2_water'};

vary_values=[[1.5,3.0,4.0,7.0,9.4,11.7,17.2]*gamma_;
            [900, 1200, 1300, 1700, 1950, 2000, 2030]*0.001;
            [90, 80, 65, 50, 42, 30, 25]*0.001; ];


%% Making simulations
clear Lobj T

selected_plot=1:7;

%Set up a structure to save simulated Z spectrum
Simulation.Sim=P;
Simulation.vary=vary;
Simulation.vary_values=vary_values;
Simulation.Zspec=zeros(size(vary_values,2), numel(P.xZspec));
Simulation.MTRasym=zeros(size(vary_values,2), numel(P.xZspec));

kplot=1;

for i=1:size(vary_values,2)
    for iv=1:numel(vary)
        v=vary{iv};
        P.(v)=vary_values(iv,i);
    end
    P=update_P(P);
    if P.analytic==1
        sol =analytic_simulation(P);
    else
        sol =numeric_simulation(P);
    end
    Simulation.Zspec(i, :)=sol.zspec;
    Simulation.MTRasym(i,:)=sol.zspec(end:-1:1)-sol.zspec;
    %Plotting
    if ismember(i,selected_plot)
        figure(1)
        yyaxis left
        T(kplot)=plot(P.xZspec, sol.zspec,'LineWidth',2.0,'LineStyle','-','Marker','none','Color',vibrant_colors(kplot)); hold on;
        %figure(2)
        yyaxis right
        posppm=(P.xZspec>=0);
        plot(P.xZspec(posppm), 100*Simulation.MTRasym(i,posppm),'LineWidth',1.5,'LineStyle','--','Marker','none','Color',vibrant_colors(kplot)); hold on;
        %xlim([0 max(P.xZspec)])
        label_list{kplot}=num2str(vary_values(1,i)/gamma_,'%.1f');
        kplot=kplot+1;
    end
end
figure(1)
Lobj=legend(T,label_list);
Lobj.Title.String='B_0 (T)';
Lobj.Location='southeast';
xlabel '\delta (ppm)'
xticks([-100,-50,-25,0,25,50,80,100])
yyaxis left
ylim([0 1])
yticks(0:0.1:1)
set(gca,'YColor','k')
ylabel 'Z'
set(gca,'XDir','reverse');
xlabel '\delta (ppm)'
yyaxis right
ylim([0 40])
yticks(0:5:15)
set(gca,'YColor','k')
ylabel 'MTR_{asym} (%)'
set(gca,'XDir','reverse');
set(gcf,'Color','w')
set(gca,'FontSize',16)
set(gca,'Box','off')

figure(3)
xexp=squeeze(vary_values(1,:)/gamma_);
yexp=100*Simulation.MTRasym(:,P.xZspec==80);
densex=[0.3:0.1:17.5];
% densey=spline(xexp,yexp,densex);
% plot(densex,densey,'LineStyle','--','Color','k','LineWidth',2.0,'Marker','none'); hold on
plot(xexp,yexp,'LineStyle','none','Color','k','LineWidth',2.0,'Marker','diamond')
ylabel 'MTR_{asym}(3 ppm) (%)'
xlabel 'B_0 (T)'
set(gcf,'Color','w')
set(gca,'Box','off')
set(gca,'FontSize',16)
ylim(100*[min(Simulation.MTRasym(:,P.xZspec==80)), max(Simulation.MTRasym(:,P.xZspec==80))])



