%% Setting up simulation scheme & generic parameters
close all

clear all

%addpath(genpath(fullfile(pwd,'qCEST_new')))

% metabolites={'glutamate_physio','glucose_physio','creatine','Phosphocreatine','Myo_inositol','Glutamine','Taurine','GABA','mAmides','NOE_-3.5ppm','MT_symm'  };
% metabolites_label={'Glu','Glc','Cr','PCr','MI','Gln','Tau','GABA','mAmides','NOE-3.5ppm','MT'};

metabolites={'glutamate_physio','creatine','mAmides','Myo_inositol','Phosphocreatine','Glutamine','Taurine','GABA','glucose_Khlebnikov'};
metabolites_label={'Glu','Cr','APT','MI','PCr','Gln','Tau','GABA','Glc'};

%more colors for plot
additional_colors=[ [0 0 .5]; [0 .5 .5]; [.5 .5 .5]; [.5 .5 0];  [0 0.5  0];[0.5 0 0]];

%% Vary one parameter
vary='B1';
vary_values=5;

xZspec=[-5:0.1:5];


Simulation.vary=vary;
Simulation.vary_values=vary_values;
Simulation.Zspec=zeros(numel(metabolites)+2,numel(vary_values), numel(xZspec));
Simulation.MTRasym=zeros(numel(metabolites)+2,numel(vary_values), numel(xZspec));
Simulation.Rex=zeros(numel(metabolites)+2,numel(vary_values), numel(xZspec));

for i=0:numel(metabolites)+2

    P=init_Sim_struct();


    if i==0
       Simulation.Sim=P;
    
    elseif i<=numel(metabolites)
            P=add_pool_model(P,metabolites{i});
    elseif i==numel(metabolites)+1
        for m=1:numel(metabolites)
            P=add_pool_model(P,metabolites{m});
        end
    elseif i==numel(metabolites)+2
        for m=1:numel(metabolites)
            P=add_pool_model(P,metabolites{m});
        end
        P=add_pool_model(P,'NOE_-3.5ppm');
    end
    
    
    P=add_pool_model(P,'MT');
    %P=add_pool_model(P,'NOE_-3.5ppm');
    P.MT_lineshape='SuperLorentzian';
    
    
    P.xZspec=xZspec;
    P.normalized=-100;    
    
    % Sets up P structure : model used for fitting
    
    %Sequence parameters
    P.FREQ          = 11.75*gamma_;       % frequency (=B0[T] * gamma)
    P.tsat=1;
    P.TR=5;
    P.Trec          = 4;                    % recovery time in s
    
    %Saturation parameters
    P.B1            = 5;                  % B1 value in uT
    P.pulsed  = 0;                    % 0 = cw saturation, 1 = pulsed saturation
    P.tp    = 1;                       % saturation time in s
    P.n     = 1;                        % choose n=1 for cw saturation
    P.shape = 'block';                  % choose 'block' for cw saturation (cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4)
    P.td=0 ;                    % choose DC=1 for cw saturation

    P.Zi=1;
    P.B1c=1.0;
    
    % 
    P.T1_water=1.95;
    P.T2_water=31.7*0.001;
    % 
    P.H_water=2*45500;
    
    %Glutamate
    P.Glu=6.5;
    P.kex_Glu=7480;
    P.T2_Glu=6.9*10^(-3);

    %Creatine
    P.Cr=4;
    P.kex_Cr=810;
    P.T2_Cr=7.1*10^(-3);

    %Myo-inositol
    P.MI=3.4;
    P.kex_MI=2090;
    P.T2_MI=22.8*10^(-3);

    %Amides
    P.APT=172;
    P.kex_APT=22;

    %PCr
    P.PCr=3.2;
    P.kex_PCr_1=67;
    P.kex_PCr_2=126;
    P.T2_PCr_1=7.8*10^(-3);
    P.T2_PCr_2=7.8*10^(-3);

    %Gln
    P.Gln=3;
    P.kex_Gln_1=17;
    P.kex_Gln_2=49;
    P.kex_Gln_3=22880;
    P.T2_Gln_1=13.8*10^(-3);
    P.T2_Gln_2=13.8*10^(-3);
    P.T2_Gln_3=13.8*10^(-3);

    %Tau
    P.Tau=12.2;
    P.kex_Tau=49600;

    %GABA
    P.GABA=1.5;
    P.kex_GABA=6900;
    P.T2_GABA=17.2*10^(-3);

    %Glc : already parametrized as Khlebnikov et al 2019

    P.MT_cutoff_ppm=2.5;
    P.dw_MT=-0.03;
    P.H_MT=(10.4/100)*P.H_water;
    P.kex_MT=13.3;
    P.H_NOEm35ppm=(4.1/100)*P.H_water;
    P.kex_NOEm35ppm=15;



    % Making simulations


    for var_val=1:numel(vary_values)
        P.(vary)=vary_values(var_val);

        P=update_P(P);
        if P.analytic==1
            sol =analytic_simulation(P);
        else
            sol =NUMERIC_SIM(P);
        end
        Simulation.Zspec(i+1,var_val, :)=sol.zspec;
        Simulation.Rex(i+1,var_val, :)=-sol.S.Rex_lorentz;
        Simulation.MTRasym(i+1,var_val,:)=sol.zspec(end:-1:1)-sol.zspec;


        %Plotting
        figure(1)
        subplot(1,2,1)
        plot(P.xZspec, sol.zspec); hold on;
        subplot(1,2,2)
        plot(P.xZspec, squeeze(Simulation.MTRasym(i+1,var_val,:))); hold on;
        xlim([0 max(P.xZspec)])
    end
    
    
    
end


%% plot contributions for all offsets


Total_CEST_signal=squeeze(Simulation.MTRasym(end,:,:));
Ref=squeeze(Simulation.Zspec(1,:,:));
%Ref=squeeze(Simulation.Zspec(end,:,end:-1:1));
Total_MTR_ref=squeeze(Simulation.Zspec(end-1,:,end:-1:1))-squeeze(Simulation.Zspec(end-1,:,:));
Real_total_MTR_asym=squeeze(Simulation.Zspec(end,:,end:-1:1))-squeeze(Simulation.Zspec(end,:,:));
%Free_MTR=squeeze(Simulation.Zspec(1,:,end:-1:1))-squeeze(Simulation.Zspec(1,:,:));

for m=1:numel(metabolites)

    MTR_metab=Ref-squeeze(Simulation.Zspec(m+1,:,:));
    % Contributions(m,:)=(Total_MTR_ref-MTR_metab);
    % MTR_ref=Ref-squeeze(Simulation.Zspec(m+1,:,:));
    % Contributions(m,:)=Total_MTR_ref-MTR_ref;
    Contributions(m,:)=MTR_metab;

end

norm=max(Total_MTR_ref);
% Contributions_norm=100*Contributions/norm;
Contributions_norm=100*Contributions;


% for ppm=1:numel(P.xZspec)
%     C=Contributions(:,ppm);
%     Contributions_norm(:,ppm)=100*C./sum(C(C>0)); %normalize by sum of positive contributions
% end

% figure(3)
% 
% tiledlayout(1,numel(vary_values))
% for i=numel(vary_values)
%    nexttile()
%    b=bar(xZspec, Contributions_norm,'stacked','FaceColor','flat');
%    title(strcat(vary,'=',num2str(vary_values(i))))
%    xlabel('\delta (ppm)')
%    ylabel('% of contribution to MTR_{asym} assesed at each ppm')
%    set(gca,'xDir','reverse')
%    xlim([0 max(xZspec)])
%    ylim([0 100])
%    
%    if numel(metabolites)>7
%         for i=8:numel(metabolites)
%             b(i).FaceColor=additional_colors(i-7,:);
%         end
%     end
% end




legend(metabolites_label,'Location','southwest')

figure(4)

tiledlayout(1,numel(vary_values))
for i=numel(vary_values)
   nexttile()
   b=bar(xZspec, Contributions_norm,'stacked','FaceColor','flat');
%    b=bar(xZspec, 100*Contributions/max(Total_MTR_ref),'stacked','FaceColor','flat');
   %title(strcat(vary,'=',num2str(vary_values(i))))
   xlabel('\delta (ppm)')
%    ylabel('% of contribution to MTR_{asym} (global normalization)')
   ylabel('Stacked MTRs (%)')
   set(gca,'xDir','reverse')
   xlim([0 max(xZspec)])
   %ylim([0 100*norm])
   
for i=1:numel(metabolites)
    if i<=4
        b(i).FaceColor=classical_4color_palette(i);
    else
       b(i).FaceColor=muted_colors(i); 
    end
end
end

index_avg=(xZspec>=1 & xZspec<=4);
for i=1:numel(metabolites)
    fprintf(' %s contribution = %.2f \n', metabolites_label{i},mean(Contributions_norm(i,index_avg)))
end


metabolites_label{numel(metabolites_label)+1}='Simulated MTR including all metabolites';
metabolites_label{numel(metabolites_label)+1}='Simulated MTR_{asym} including all metabolites + NOE^1';


hold on
plot(xZspec(xZspec>=0),100*Total_MTR_ref(xZspec>=0),'k','LineStyle','--','LineWidth',1.5)

plot(xZspec(xZspec>=0),100*Real_total_MTR_asym(xZspec>=0),'r','LineStyle','--','LineWidth',1.5)

legend(metabolites_label,'Location','northwest')

set(gca,'Box','off')
set(gca,'FontSize',14)
set(gcf,'Color','w')

title 'Simulation of contributions of metabolites to MTR'

%% plot contributions at one offset

offset=3.0; %choose offset

Total_MTR_approx=sum(Contributions,1);
Normed_contrib=100*Contributions./Total_MTR_approx;
Contributions_offset=Normed_contrib(:,find_nearest(P.xZspec, offset));
%Metabolite contribution is (Total_simulation-Simulation_without_metabolite)
%Contributions_offset=(Total_CEST_signal_offset-(Ref-squeeze(Simulation.Zspec(2:end-1,:,find_nearest(P.xZspec, offset)))));
%Contributions_offset=100*Contributions_offset/sum(Contributions_offset(Contributions_offset>0)); %normalize by sum of positive contributions
% Contributions_offset=(Total_MTR_ref_offset-(Ref-squeeze(Simulation.Zspec(2:end-1,:,find_nearest(P.xZspec, offset)))));
% Contributions_offset=100*Contributions_offset/Total_MTR_ref_offset;

figure(2)
% subplot(1,2,1)
b=bar(vary_values, Contributions_offset, 'stacked','FaceColor','flat');
xlabel('Stacked MTRs')
ylabel('Contributions (% of summed MTR)')

for i=1:numel(metabolites)
    if i<=4
        b(i).FaceColor=classical_4color_palette(i);
    else
       b(i).FaceColor=muted_colors(i); 
    end
end

title(strcat('Contributions at \delta =',num2str(offset),' ppm'));

% subplot(1,2,2)
% b=bar(vary_values, Contributions_offset, 'stacked','FaceColor','flat');
% title('Only positive contributions')
% xlabel(vary)
% ylabel('% of contributions to MTRasym')
% % ylim([0 100])

% if numel(metabolites)>7
%     for i=8:numel(metabolites)
%         b(i).FaceColor=additional_colors(i-7,:);
%     end
% end

legend(metabolites_label)


for i=1:numel(metabolites)
    fprintf(' %s contribution = %.2f \n', metabolites_label{i},mean(Contributions_offset(i,:)))
end

xticks([])
set(gca,'Box','off')
set(gca,'FontSize',14)
set(gcf,'Color','w')


%% plot contributions for all offsets


% Total_CEST_signal=squeeze(Simulation.MTRasym(end,:,:));
% Ref=squeeze(Simulation.Zspec(1,:,:));
% %Ref=squeeze(Simulation.Zspec(end,:,end:-1:1));
Total_Rex_ref=squeeze(Simulation.Rex(end-1,:,:));
% Real_total_MTR_asym=squeeze(Simulation.Zspec(end,:,end:-1:1))-squeeze(Simulation.Zspec(end,:,:));
% %Free_MTR=squeeze(Simulation.Zspec(1,:,end:-1:1))-squeeze(Simulation.Zspec(1,:,:));

for m=1:numel(metabolites)

%     MTR_metab=Ref-squeeze(Simulation.Zspec(m+1,:,:));
    % Contributions(m,:)=(Total_MTR_ref-MTR_metab);
    % MTR_ref=Ref-squeeze(Simulation.Zspec(m+1,:,:));
    % Contributions(m,:)=Total_MTR_ref-MTR_ref;
    Contributions_Rex(m,:)=Simulation.Rex(m+1,:);

end

norm_Rex=1;%max(Total_MTR_ref);
Contributions_norm_Rex=Contributions_Rex/norm_Rex;
%Contributions_norm=100*Contributions_Rex;





legend(metabolites_label,'Location','southwest')

figure(5)

tiledlayout(1,numel(vary_values))
for i=numel(vary_values)
   nexttile()
   b=bar(xZspec, Contributions_norm_Rex,'stacked','FaceColor','flat');
%    b=bar(xZspec, 100*Contributions/max(Total_MTR_ref),'stacked','FaceColor','flat');
   %title(strcat(vary,'=',num2str(vary_values(i))))
   xlabel('\delta (ppm)')
%    ylabel('% of contribution to MTR_{asym} (global normalization)')
   ylabel('Stacked R_{ex}^j (s^{-1})')
   set(gca,'xDir','reverse')
   
    xlim([0 max(xZspec)])
   %ylim([0 100*norm])
   
for i=1:numel(metabolites)
    if i<=4
        b(i).FaceColor=classical_4color_palette(i);
    else
       b(i).FaceColor=muted_colors(i); 
    end
end
end

% index_avg=(xZspec>=1 & xZspec<=4);
% for i=1:numel(metabolites)
%     fprintf(' %s contribution = %.2f \n', metabolites_label{i},mean(Contributions_norm(i,index_avg)))
% end


metabolites_label{numel(metabolites_label)+1}='Simulated R_{ex}^{total} including all metabolites';
% metabolites_label{numel(metabolites_label)+1}='Simulated MTR_{asym} including all metabolites + NOE^1';
% 
% 
hold on
plot(xZspec(xZspec>=0),Total_Rex_ref(xZspec>=0),'k','LineStyle','--','LineWidth',1.5)
% % 
% % plot(xZspec(xZspec>=0),100*Real_total_MTR_asym(xZspec>=0),'r','LineStyle','--','LineWidth',1.5)

legend(metabolites_label,'Location','northwest')

set(gca,'Box','off')
set(gca,'FontSize',14)
set(gcf,'Color','w')

title 'Simulation of contributions of metabolites to R_{ex}^{total}'

offset=3.0; %choose offset

Total_Rex=sum(Contributions_Rex,1);
Normed_contrib_Rex=100*Contributions_Rex./Total_Rex;
Contributions_offset_Rex=Normed_contrib_Rex(:,find_nearest(P.xZspec, offset));
%Metabolite contribution is (Total_simulation-Simulation_without_metabolite)
%Contributions_offset=(Total_CEST_signal_offset-(Ref-squeeze(Simulation.Zspec(2:end-1,:,find_nearest(P.xZspec, offset)))));
%Contributions_offset=100*Contributions_offset/sum(Contributions_offset(Contributions_offset>0)); %normalize by sum of positive contributions
% Contributions_offset=(Total_MTR_ref_offset-(Ref-squeeze(Simulation.Zspec(2:end-1,:,find_nearest(P.xZspec, offset)))));
% Contributions_offset=100*Contributions_offset/Total_MTR_ref_offset;

figure(6)
% subplot(1,2,1)
b=bar(vary_values, Contributions_offset_Rex, 'stacked','FaceColor','flat');
xlabel('Stacked R_{ex}^j')
ylabel('Contributions (% of R_{ex}^{total})')

for i=1:numel(metabolites)
    if i<=4
        b(i).FaceColor=classical_4color_palette(i);
    else
       b(i).FaceColor=muted_colors(i); 
    end
end

title(strcat('Contributions to R_{ex}^{total} at \delta =',num2str(offset),' ppm'));

% subplot(1,2,2)
% b=bar(vary_values, Contributions_offset, 'stacked','FaceColor','flat');
% title('Only positive contributions')
% xlabel(vary)
% ylabel('% of contributions to MTRasym')
% % ylim([0 100])

% if numel(metabolites)>7
%     for i=8:numel(metabolites)
%         b(i).FaceColor=additional_colors(i-7,:);
%     end
% end

legend(metabolites_label)


% for i=1:numel(metabolites)
%     fprintf(' %s contribution = %.2f \n', metabolites_label{i},mean(Contributions_offset(i,:)))
% end

xticks([])
set(gca,'Box','off')
set(gca,'FontSize',14)
set(gcf,'Color','w')
