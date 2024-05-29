function plot_fit_results_table(Results)
%Plot into a table all the parameters of the simulation

%Table for general parameters results
header={'T1 (s)', 'T2 (s)', 'Proton fraction f', 'Exchange rate /w water (Hz)' , 'Chemical shift (ppm)'};
columnT1= [P.T1_water; 1/Sim.R1B; 1/Sim.R1C];
columnT2= [1/Sim.R2A; 1/Sim.R2B; 1/Sim.R2C];
columnfi= [ 1.0; Sim.fB; Sim.fC];
columnkiA= [ NaN; Sim.kBA; Sim.kCA];
columndwi= [Sim.dwA; Sim.dwB; Sim.dwC];

Tpools=table(columnT1, columnT2, columnfi, columnkiA, columndwi, 'RowNames', Sim.pool_names);
Tpools.Properties.VariableNames=header;

%Table for sequence parameters
header={'B1 (uT)', 'Tsat (s)', 'Pulsed (=0) or continous (=1) saturation?' };

Tseq=table(Sim.B1, Sim.tp, Sim.pulsed);
%Tseq=table(Sim.B1, Sim.tp);
Tseq.Properties.VariableNames=header;   


f=figure(445);
uitable(f,'Data',Tpools{:,:},'ColumnName',Tpools.Properties.VariableNames,...
    'RowName',Tpools.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

uitable(f,'Data',Tseq{:,:},'ColumnName',Tseq.Properties.VariableNames,...
   'RowName','Acquisition','Units', 'Normalized', 'Position',[0, -0.3, 1, 1]);


end