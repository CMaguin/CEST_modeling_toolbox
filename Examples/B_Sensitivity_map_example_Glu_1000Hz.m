%% Setting up simulation scheme & generic parameters
close all

clear all
P=init_Sim_struct();




%Sequence parameters
P.xZspec=[-6:0.05:6];

P.FREQ=11.75*gamma_; %B0=11.7T
P.Zi=1.0;
P.B1=5;		 % the saturation B1 in µT

P=add_pool_model(P,'glutamate_physio');

P.T1_water=2.0;
P.T2_water=0.03;

P.Glu=15;
P.kex_Glu=1000;
P.dw_Glu=3.0;
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
vary={'B1','tp'};

vary_values_B1=0:0.1:10;
vary_values_tsat=0:0.1:10;
vary_values=zeros(numel(vary_values_B1),numel(vary_values_tsat),2);
for i=1:numel(vary_values_B1)
    for j=1:numel(vary_values_tsat)
        vary_values(i,j,:)=[vary_values_B1(i),vary_values_tsat(j)];
    end
end

sv=size(vary_values);

sv=size(vary_values);
%% Making simulations

%Set up a structure to save simulated Z spectrum
Simulation.Sim=P;
Simulation.vary=vary;
Simulation.vary_values=vary_values;
Simulation.Zspec=zeros(sv(1)*sv(2), numel(P.xZspec));
Simulation.MTRasym=zeros(sv(1)*sv(2), numel(P.xZspec));

kplot=1;

tic

selected_plotB1=[1,3,5];
selected_plottp=[1,4];

for i=1:sv(1)*sv(2)
    
    
    for vi=1:numel(vary)
        [v(1),v(2)]=ind2sub([numel(vary_values_B1),numel(vary_values_tsat)],i);
        varyvariable=vary{vi};
        P.(varyvariable)=vary_values(v(1),v(2),vi);
        fprintf(' %s = %1f \n',varyvariable,vary_values(v(1),v(2),vi))
        cB1=vary_values(v(1),v(2),1);
        ctp=vary_values(v(1),v(2),2);
    end
    P=update_P(P);
    if P.analytic==1
        sol =analytic_simulation(P);
    else
        sol =numeric_simulation(P);
    end
    Simulation.Zspec(i, :)=sol.zspec;
    Simulation.MTRasym(i,:)=sol.zspec(end:-1:1)-sol.zspec;
    
    
    if ismember(cB1,selected_plotB1)
        if ismember(ctp,selected_plottp)
        figure(1)
        yyaxis left
        T(kplot)=plot(P.xZspec, sol.zspec,'LineWidth',2.0,'LineStyle','-','Marker','none','Color',diverging_colors_2(kplot,10)); hold on;
        %figure(2)
        yyaxis right
        posppm=(P.xZspec>=0);
        plot(P.xZspec(posppm), 100*Simulation.MTRasym(i,posppm),'LineWidth',1.5,'LineStyle','--','Marker','none','Color',diverging_colors_2(kplot,10)); hold on;
        %xlim([0 max(P.xZspec)])
        label_list{kplot}=strcat('B_1=',num2str(cB1,'%.0f') ,'µT; t_{sat}=',num2str(ctp,'%.0f'),'s');
        kplot=kplot+1;
        end
    end
end

toc

Lobj=legend(T,label_list);
Lobj.Location='southeast';
yyaxis left
ylim([0 1])
yticks(0:0.1:1)
set(gca,'YColor','k')
ylabel 'Z'
set(gca,'XDir','reverse');
xlabel '\delta (ppm)'
yyaxis right
xlabel '\delta (ppm)'
ylim([0 40])
yticks(0:5:15)
set(gca,'YColor','k')
ylabel 'MTR_{asym} (%)'
set(gca,'XDir','reverse');
set(gcf,'Color','w')
set(gca,'FontSize',16)
set(gca,'Box','off')
title 'Examples of Z-spectra'


%% Genrating sensitivity map
f=figure(12);


CEST_signal=100*Simulation.MTRasym(:,P.xZspec==P.dw_Glu);
CEST_signal=reshape(CEST_signal, sv(1),sv(2));

imagesc(vary_values_tsat,vary_values_B1,CEST_signal); hold on
c=colorbar();
%clim([0 15])
colormap(cmocean('haline'))
c.Label.String = strcat('MTR_{asym} @', num2str(P.dw_Glu),' ppm (%)');
c.Limits=[0 12];
% xlabel(vary{2})
% ylabel(vary{1})
xlabel('t_{sat} (s)')
ylabel('B_1 (µT)')
set(gca,'YDir','normal')
xlim([0 6])
yticks(0:10)
xticks(0:10)


for j=1:numel(vary_values_tsat)
    [maxMTR,iB1opt]=max(CEST_signal(:,j));
    B1opt(j)=vary_values_B1(iB1opt);
end
plot(vary_values_tsat,B1opt,'LineStyle','--','Color',[0.25 0.25 0.25],'LineWidth',2.0);
legend('B_1^{max}')

set(gcf,'Color','w')
set(gca,'FontSize',16)
set(gca,'Box','off')

