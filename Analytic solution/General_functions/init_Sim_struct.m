function P=init_Sim_struct()
    P.n_cest_pools=0;
    P.pool_names={};
    P.H=[];
    P.kex=[];
    P.dw=[];
    P.T1=[];
    P.T2=[];

    %By default no MT and analytic Rex solution
    P.MT=0;
    P.analytic=1;
    P.Rex_sol       = 'Lorentz'; 

    %Init water pool with default values
    P.model={'Water'};
    P.dw_water=0;
    P.H_water=2*55.6*1000; %water proton concentration in mM
    P.T1_water=2.2; %T1 in s
    P.T2_water=0.1; %T2 in s
    
    
    P.play_readout=0; %by default readout sequence is not taken into account in the simulation
    
    P.WaterResCorr=0;
end