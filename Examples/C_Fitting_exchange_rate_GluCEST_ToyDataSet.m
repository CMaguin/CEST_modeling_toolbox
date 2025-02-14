clear all
clc
close all

%addpath(genpath(pwd))
%% load data
ToyDataSetpath=fullfile(pwd, 'ToyDataSet');

B1values=[1:10];

Glu=[100.1,100.1,99.89,99.89,99.93,99.93]; %Exact experimental concentration of [Glu]
Buffer={'H2O_22C','H2O_37C','HEPES_22C','HEPES_37C','PBS_22C','PBS_37C'};
Ntubes=numel(Buffer);


xZspec_full=load(fullfile(ToyDataSetpath,'H2O_22C','freq_ppm_1.mat')).freq_ppm;
xind_red=(abs(xZspec_full)<=5);
xZspec=xZspec_full(xind_red);
if isrow(xZspec)
xZspec = xZspec';
end
xxx=xZspec;

load(fullfile(ToyDataSetpath,'T1_study.mat'))
load(fullfile(ToyDataSetpath,'T2_study.mat'))


noffsets=numel(xZspec);
Zspec=zeros(numel(Buffer),numel(B1values), noffsets);

xZspec_interp=-5:0.1:5;
MTRasym_avg=zeros(numel(Buffer),numel(B1values), numel(xZspec_interp));

for b=1:numel(Buffer)
    Buffer_name=Buffer{b};
    datapath=fullfile(ToyDataSetpath,Buffer_name);
    for ZNum=1:numel(B1values)
        
        Data=load(strcat(datapath, '\Zspectrum_',num2str(ZNum),'.mat'));
        f=fieldnames(Data);
        Data=Data.(f{1});

        index_M0=3;
        ppmM0=xZspec_full(index_M0);
    
        Data_norm=Data./Data(index_M0);

       
    
        Zspec(b,ZNum,:)=Data_norm(xind_red);

        
        Zspec_interp(b,ZNum,:)=interp1(xZspec,squeeze(Zspec(b,ZNum,:)),xZspec_interp,'pchip');

        MTRasym_avg(b,ZNum,:)=Zspec_interp(b,ZNum,end:-1:1)-Zspec_interp(b,ZNum,:);

    end
end


%% Plot data for one phantom at multiple B1 values

clear T labels

selected_B1=1:10;
selected_tubes=1:6;

for i=selected_tubes
    figure(i)
    for ZNum=selected_B1
        yyaxis left 
        T(ZNum)=plot(xZspec_interp,squeeze(Zspec_interp(i,ZNum,:)),'LineStyle','-','Color',diverging_colors_2(ZNum,numel(selected_B1)),'Marker','none','LineWidth',1.5);
        hold on
        yyaxis right
        pos_ppm=(xZspec_interp>=0.2);
        plot(xZspec_interp(pos_ppm),100*squeeze(MTRasym_avg(i,ZNum,pos_ppm)),'LineStyle','--','Color',diverging_colors_2(ZNum,numel(selected_B1)),'Marker','none','LineWidth',1.5);
        hold on;
        labels{ZNum}=strcat(num2str(B1values(ZNum)), ' µT');
    end


set(gca,'XDir','reverse')
set(gcf,'Color','w')
xlabel('\delta (ppm)')
yyaxis left
ylabel 'Z'
ylim([0 1])
yticks(0:0.2:1)
set(gca,'YColor','k')
yyaxis right
ylabel 'MTR_{asym} (%)'
ylim([0 2*max(100*squeeze(MTRasym_avg(i,8,pos_ppm)))])
yticks(100*[0:0.05:0.8])
set(gca,'YColor','k')
set(gca,'Box','off')
set(gcf,'Color','w')

legend(T,labels,'Location','southeast','Interpreter','none','FontSize',12);

title(Buffer{i},'Interpreter','none')

end

%% Plot Z-spectra of every buffer for each temperature at 5 µT 
clear T labels

selected_B1=5;

%at Room Temperature
figure(8)
title 'Room Temperature ([Glu] = 100 mM)'
selected_tubes=[1,3,5];
kc=0;
for i=selected_tubes
    kc=kc+1;
    for ZNum=selected_B1
        yyaxis left
        T(kc)=plot(xZspec_interp,squeeze(Zspec_interp(i,ZNum,:)),'LineStyle','-','Color',classical_4color_palette(kc),'Marker','none','LineWidth',2.0);
        hold on
        yyaxis right
        pos_ppm=(xZspec_interp>=0.2);
        plot(xZspec_interp(pos_ppm),100*squeeze(MTRasym_avg(i,ZNum,pos_ppm)),'LineStyle','-','Color',classical_4color_palette(kc),'Marker','none','LineWidth',1.0);
        hold on;
        labels{kc}=strcat('100 mM Glu in  ',Buffer{i});
    end

    %labels={'100 mM Glu in water','100 mM Glu in HEPES','100 mM Glu in PBS'};

set(gca,'XDir','reverse')
set(gcf,'Color','w')
xlabel('\delta (ppm)')
yyaxis left
ylabel 'Z'
ylim([0 1.02])
yticks(0:0.2:1)
set(gca,'YColor','k')
yyaxis right
ylabel 'MTR_{asym} (%)'
ylim([0 1.7]*100)
yticks(100*[0:0.1:0.6])
set(gca,'YColor','k')
set(gca,'Box','off')
set(gcf,'Color','w')

legend(T,labels,'Location','southeast','Interpreter','none','FontSize',12);

end

%at Room Temperature
figure(9)
title '37°C ([Glu] = 100 mM)'
selected_tubes=[2,4,6];
kc=0;
for i=selected_tubes
    kc=kc+1;
    for ZNum=selected_B1
        yyaxis left
        T2(kc)=plot(xZspec_interp,squeeze(Zspec_interp(i,ZNum,:)),'LineStyle','-','Color',classical_4color_palette(kc),'Marker','none','LineWidth',2.0);
        hold on
        yyaxis right
        pos_ppm=(xZspec_interp>=0.2);
        plot(xZspec_interp(pos_ppm),100*squeeze(MTRasym_avg(i,ZNum,pos_ppm)),'LineStyle','-','Color',classical_4color_palette(kc),'Marker','none','LineWidth',1.0);
        hold on;
        labels2{kc}=strcat('100 mM Glu in  ',Buffer{i});
    end

    %labels={'100 mM Glu in water','100 mM Glu in HEPES','100 mM Glu in PBS'};

set(gca,'XDir','reverse')
set(gcf,'Color','w')
xlabel('\delta (ppm)')
yyaxis left
ylabel 'Z'
ylim([0 1.02])
yticks(0:0.2:1)
set(gca,'YColor','k')
yyaxis right
ylabel 'MTR_{asym} (%)'
ylim([0 1.7]*100)
yticks(100*[0:0.1:0.6])
set(gca,'YColor','k')
set(gca,'Box','off')
set(gcf,'Color','w')

legend(T2,labels2,'Location','southeast','Interpreter','none','FontSize',12);

end



%% kex fit
selected_exp=1:Ntubes;
selected_B1=[1:8]; % removed some B1 values
for selected_tubes=selected_exp

    Z=squeeze(Zspec(selected_tubes,selected_B1,:));
    P=init_Sim_struct();


    P.xZspec=xZspec;
    P.normalized=-20;
    
    P.Glu=Glu(selected_tubes);
    P=add_pool_model(P,'glutamate_physio');



    % Sets up P structure : model used for fitting
    
    %Sequence parameters
    P.FREQ          = 11.75*gamma_;       % frequency (=B0[T] * gamma)
    P.tsat=1;
    P.B1            = B1values(1);                  % B1 value in ÂµT
    P.TR=5;
    P.Trec          = 4;                    % recovery time in s
    P.Zi            = 1.0;                    % initial magnetisation 
    P.M0corr        = 1.0;                    % correction on M0
    
    % P.Offsets_order='successive'; %if you want to take into account incomplete relaxation between two offsets 

    P.pulsed  = 0;                    % 0 = cw saturation, 1 = pulsed saturation


    P.tp    = 1;                       % saturation time in s
    P.n     = 1;                        % choose n=1 for cw saturation
    P.shape = 'block';                  % choose 'block' for cw saturation (cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4)
    P.td=0 ;                    % choose DC=1 for cw saturation


    P.T1_water=Inferred_T1_H2O(selected_tubes)/1000;%T1 estimation (in s)
    P.T2_water=Inferred_T2_H2O(selected_tubes)/1000;%T2 estimation (in s)
    
    P=update_P(P);

    P.WaterResCorr=1;

    % Sets up OPTIM structure : which variables are allowed to vary?
    
    OPTIM=init_OPTIM_struc(1, numel(B1values(selected_B1)));
    
    OPTIM.fit_type=0; % fit only Z-spectrum

    %Set up varying parameters
    OPTIM.vary_indiv={'B1'};
    OPTIM.varyval_indiv=B1values(selected_B1);


    % OPTIM=add_voxel_var(OPTIM, 'B1c', '1', '0.85', '1.15'); % if you want B1 correction as a free variable
    % OPTIM=add_indiv_var(OPTIM, 'M0corr','1.0','0.95','1.05'); % if you want M0 correction as a free variable

    %OPTIM=add_voxel_var(OPTIM, 'Glu', 'P.Glu', 'P.Glu*0.95', 'P.Glu*1.05'); % if you want to add freedom on [Glu] 
    OPTIM=add_voxel_var(OPTIM, 'kex_Glu', '7000', '0', '2000000');
    OPTIM=add_voxel_var(OPTIM, 'T1_water', 'P.T1_water', '1.8', '5');
    OPTIM=add_voxel_var(OPTIM, 'T2_water', 'P.T2_water', '0.05', '0.5');

    %OPTIM=add_voxel_var(OPTIM, 'dw_water', 'P.dw_water','P.dw_water-0.1','P.dw_water+0.1'); %if you want to add freedom on B0 correction

    
    OPTIM=define_start_and_bounds(OPTIM,P);



    % Call fitting function 
    options=optimset('MaxIter',5000);
    tic
    [Results,Fit, newP]=fit_data(P,Z,OPTIM,options);
    toc
    
    Resulting_fit(selected_tubes)=Fit;
    R2(selected_tubes)=Results.R_sq;
    Inferred_res{selected_tubes}=Results;
    Inferred_kex(selected_tubes)=Results.voxel_val(:,strcmp(Results.voxel_vars,'kex_Glu'));
    Error_on_kex(selected_tubes)=Results.ci_voxel(:,strcmp(Results.voxel_vars,'kex_Glu'),1);
    fprintf('Phantom made in %s \n',Buffer{selected_tubes})
    fprintf('R^2 = %.5f \n',Fit.Rsq);
    fprintf('Estimated kex_Glu = %.4f +- %.4f Hz \n',Inferred_kex(selected_tubes),Results.ci_voxel(:,strcmp(Results.voxel_vars,'kex_Glu'),1))
    fprintf('Estimated T1w = %.2f +- %.2f s \n',Results.voxel_val(:,strcmp(Results.voxel_vars,'T1_water')),Results.ci_voxel(:,strcmp(Results.voxel_vars,'T1_water'),1))
    fprintf('Estimated T2w = %.4f +- %.4f s \n',Results.voxel_val(:,strcmp(Results.voxel_vars,'T2_water')),Results.ci_voxel(:,strcmp(Results.voxel_vars,'T2_water'),1))
    %fprintf('Estimated B0 corr = %.4f +- %.4f ppm \n',Results.voxel_val(:,strcmp(Results.voxel_vars,'dw_water')),Results.ci_voxel(:,strcmp(Results.voxel_vars,'dw_water'),1))
    %fprintf('Estimated B1 corr = %.4f +- %.4f ppm \n',Results.voxel_val(:,strcmp(Results.voxel_vars,'B1c')),Results.ci_voxel(:,strcmp(Results.voxel_vars,'B1c'),1))
    

    
    %plotting fit
    figure(136+selected_tubes); % i like the number 136
    iB1=1;
    for ZNum=selected_B1
        plot(xZspec,Z(iB1,:),'Color',diverging_colors_2(ZNum,10),'LineStyle','none','Marker','+','MarkerSize',10); hold on    
        Tx(iB1)=plot(Fit.densex,squeeze(Fit.fitdense(:,iB1,:)),'Color',diverging_colors_2(ZNum,10),'LineStyle','-','Marker','none','LineWidth',1.5); hold on
        lab{iB1}=strcat(num2str(B1values(ZNum)),' µT');
        iB1=iB1+1;
        
    end
    legend(Tx,lab,'Location','southeast');
    set(gca,'XDir','reverse')
    set(gcf,'Color','w')
    xlabel('\delta (ppm)')
    ylabel 'Z'
    title(strcat('Fit of multi-B1 CEST data of a 100 mM Glu phantom in ',Buffer{selected_tubes}),'Interpreter','none')
end