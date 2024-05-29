function P=add_pool_model(P, metabolite_name,option)

P.model{numel(P.model)+1}=metabolite_name;
switch metabolite_name

    case 'MT'  
        P.MT=1;
        P.MT_lineshape='SuperLorentzian'; %default
        P.H_MT=0.05*2*55.6*1000;    % 5%
        P.kex_MT=40;      %Hz   
        %P.dw_MT=0; %symmetric case
        %P.dw_MT=-2.6;             %asymmetric MT according to Zaiss (?)
        P.dw_MT=-2.34;             %asymmetric MT according to Hua 2007
        P.T1_MT=1; 
        P.T2_MT=1e-5; %10Âµs    
        
    case 'MT_symm'  
        P.MT=1;
        P.MT_lineshape='SuperLorentzian'; %default
        P.H_MT=0.1*2*55.6*1000;    %10%
        P.kex_MT=40;      %Hz   
        P.dw_MT=0;             
        P.T1_MT=1; 
        P.T2_MT=1e-5; %10Âµs  

    case 'MT_flex'  
        P.MT=1;
        P.MT_lineshape='Flexible'; %default
        P.H_MT=0.1*2*55.6*1000;    %10%
        P.kex_MT=40;      %Hz   
        P.dw_MT=0;             
        P.T1_MT=1; 
        P.T2_MT=1e-5; %10Âµs  
        P.MT_flex_fct=@(dMT) 0*dMT; %fill this up with our experimentally found flex lineshape
        
    case 'NOE_-3.5ppm' %Zhang 2017
        P.pool_names{1+P.n_cest_pools}='NOEm35ppm';
        P.H=[P.H                778];    %fs=0.007
        P.kex=[P.kex            50];      %Hz 
        P.dw=[P.dw              -3.5];   %ppm
        P.T1=[P.T1      1.5]; 
        P.T2=[P.T2      0.0005]; 
        P.n_cest_pools=P.n_cest_pools+1;
        
    case 'NOE_-1.6ppm' %Zhang 2017
        P.pool_names{1+P.n_cest_pools}='NOEm16ppm';
        P.H=[P.H                333];    %fs=0.003
        P.kex=[P.kex            50];      %Hz 
        P.dw=[P.dw              -1.6];   %ppm
        P.T1=[P.T1      1.5]; 
        P.T2=[P.T2      0.001]; 
        P.n_cest_pools=P.n_cest_pools+1;
        
    case 'NOE_-1.75ppm' %Van Zilj 2018
        P.pool_names{1+P.n_cest_pools}='NOEm175ppm';
        P.H=[P.H                100];    
        P.kex=[P.kex            16];      %Hz 
        P.dw=[P.dw              -1.75];   %ppm
        P.T1=[P.T1      P.T1_water]; 
        P.T2=[P.T2      0.005]; 
        P.n_cest_pools=P.n_cest_pools+1;
        
    case 'NOE_-2.25ppm' %Van Zilj 2018
        P.pool_names{1+P.n_cest_pools}='NOEm225ppm';
        P.H=[P.H                100];   
        P.kex=[P.kex            16];      %Hz 
        P.dw=[P.dw              -2.25];   %ppm
        P.T1=[P.T1      P.T1_water]; 
        P.T2=[P.T2      0.005]; 
        P.n_cest_pools=P.n_cest_pools+1;
        
    case 'NOE_-2.75ppm' %Van Zilj 2018
        P.pool_names{1+P.n_cest_pools}='NOEm275ppm';
        P.H=[P.H                100];   
        P.kex=[P.kex            16];      %Hz 
        P.dw=[P.dw              -2.75];   %ppm
        P.T1=[P.T1      P.T1_water]; 
        P.T2=[P.T2      0.005]; 
        P.n_cest_pools=P.n_cest_pools+1;
        
    case 'NOE_-3.25ppm' %Van Zilj 2018
        P.pool_names{1+P.n_cest_pools}='NOEm325ppm';
        P.H=[P.H                100];   
        P.kex=[P.kex            16];      %Hz 
        P.dw=[P.dw              -3.25];   %ppm
        P.T1=[P.T1      P.T1_water]; 
        P.T2=[P.T2      0.005]; 
        P.n_cest_pools=P.n_cest_pools+1;
        
    case 'NOE_-3.75ppm' %Van Zilj 2018
        P.pool_names{1+P.n_cest_pools}='NOEm375ppm';
        P.H=[P.H                100];   
        P.kex=[P.kex            16];      %Hz 
        P.dw=[P.dw              -3.75];   %ppm
        P.T1=[P.T1      P.T1_water]; 
        P.T2=[P.T2      0.005]; 
        P.n_cest_pools=P.n_cest_pools+1;
        
    case 'mAmides'  %as given in Van Zilj 2019
        if ~isfield(P,'APT'); P.APT=166; end    %Zhang 2017
        P.pool_names{1+P.n_cest_pools}='APT';
        P.H=[P.H                P.APT];    %2 exchangeable protons in amide
        P.kex=[P.kex            22];      %Hz Khlebnikov 2018
        P.dw=[P.dw              3.5];         %ppm
        P.T1=[P.T1      P.T1_water]; %by default T1 is the same as water
        P.T2=[P.T2      0.002]; %Zhang 2017
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=1*P.APT;'));
        
        P.n_cest_pools=P.n_cest_pools+1;
    
    case 'mAmides_secondary'  %copy, for another amide exchanging pool
        if ~isfield(P,'APT_2'); P.APT_2=100; end    
        P.pool_names{1+P.n_cest_pools}='mAmides_2';
        P.H=[P.H                P.APT_2];    %2 exchangeable protons in amide
        P.kex=[P.kex            150];      %Hz Khlebnikov 2018
        P.dw=[P.dw              3.5];         %ppm
        P.T1=[P.T1      P.T1_water]; %by default T1 is the same as water
        P.T2=[P.T2      0.002]; %Zhang 2017
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=1*P.APT_2;'));
        
        P.n_cest_pools=P.n_cest_pools+1; 
        
    case 'Amines' %a global pool of amines resonating at 1.9ppm
        P.pool_names{1+P.n_cest_pools}='Amines';
        P.H=[P.H                10];    
        P.kex=[P.kex            1000];      %Hz 
        P.dw=[P.dw              1.9];         %ppm
        P.T1=[P.T1      P.T1_water]; %by default T1 is the same as water
        P.T2=[P.T2      0.0078];         
        
        P.n_cest_pools=P.n_cest_pools+1; 
        
    case 'Amines_qp' %as proposed in Quantiphyse CEST module
        P.pool_names{1+P.n_cest_pools}='Amines';
        P.H=[P.H                10];    
        P.kex=[P.kex            500];      %Hz 
        P.dw=[P.dw              2.8];         %ppm
        P.T1=[P.T1      P.T1_water]; %by default T1 is the same as water
        P.T2=[P.T2      0.00025];         
        
        P.n_cest_pools=P.n_cest_pools+1; 
        
    case 'creatine'  %as given in Van Zilj 2019
        if ~isfield(P,'Cr'); P.Cr=6; end    %if not given, init [Cr] at 6mM
        P.pool_names{1+P.n_cest_pools}='Cr';
        P.H=[P.H                4*P.Cr];    %4 exchangeable protons at pH=7
        %P.kex=[P.kex            1100];      %Hz Van Zilj  
        P.kex=[P.kex            810];      %Hz Khlebnikov 2018
%         P.kex=[P.kex            713];      %Hz Wermter 32°C
        P.dw=[P.dw              2];         %ppm
        P.T1=[P.T1      P.T1_water]; %by default T1 is the same as water
        P.T2=[P.T2      0.0078]; %Khlebnikov 2018
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=4*P.Cr;'));
        P.n_cest_pools=P.n_cest_pools+1;
        

        
    case 'Myo_inositol'
        if ~isfield(P,'MI'); P.MI=5; end    %Khlebnikov 2018
        P.pool_names{1+P.n_cest_pools}='MI';
        P.H=[P.H                6*P.MI];    %6 exchangeable protons at pH=7
        P.kex=[P.kex            2090];      %Hz Khlebnikov 2018 at pH=7
        P.dw=[P.dw              1.0];         %ppm
        P.T1=[P.T1      P.T1_water]; %by default T1 is the same as water
        P.T2=[P.T2      0.023]; %Khlebnikov 2018
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=6*P.MI;'));
        
        P.n_cest_pools=P.n_cest_pools+1; 
        
    case 'GABA'
        if ~isfield(P,'GABA'); P.GABA=1.5; end    %Khlebnikov 2018
        P.pool_names{1+P.n_cest_pools}='GABA';
        P.H=[P.H                3*P.GABA];    %3 exchangeable protons at pH=7
        P.kex=[P.kex            6900];      %Hz Khlebnikov 2018 at pH=7
        P.dw=[P.dw              2.91];         %ppm
        P.T1=[P.T1      P.T1_water]; %by default T1 is the same as water
        P.T2=[P.T2      0.017]; %Khlebnikov 2018
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=3*P.GABA;'));
        P.n_cest_pools=P.n_cest_pools+1; 
        
    case 'Taurine'
        if ~isfield(P,'Tau'); P.Tau=1.5; end    %Khlebnikov 2018
        P.pool_names{1+P.n_cest_pools}='Tau';
        P.H=[P.H                3*P.Tau];    %3 exchangeable protons at pH=7
        P.kex=[P.kex            49600];      %Hz Khlebnikov 2018 at pH=7
        P.dw=[P.dw              3.18];         %ppm
        P.T1=[P.T1      P.T1_water]; %by default T1 is the same as water
        P.T2=[P.T2      0.01]; %arbitrary (but doesnt really matte for this value of kex)
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=3*P.Tau;'));
        P.n_cest_pools=P.n_cest_pools+1; 
        
    case 'Phosphocreatine'
        if ~isfield(P,'PCr'); P.PCr=6; end    %Khlebnikov 2018
        P.pool_names{1+P.n_cest_pools}='PCr_1';
        P.pool_names{2+P.n_cest_pools}='PCr_2';
        P.H=[P.H                2*P.PCr      1*P.PCr];    %form at pH=7
        P.kex=[P.kex            67      128];      %Hz Khlebnikov 2018 at pH=7
        P.dw=[P.dw              1.93    2.64];         %ppm
        P.T1=[P.T1      P.T1_water  P.T1_water]; %by default T1 is the same as water
        P.T2=[P.T2      0.008       0.008]; %Khlebnikov 2018
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=2*P.PCr;'));
        P=add_constrained_relation(P,strcat('P.H(',num2str(2+P.n_cest_pools),')=1*P.PCr;'));
        P.n_cest_pools=P.n_cest_pools+2; 
        
    case 'Glutamine'
        if ~isfield(P,'Gln'); P.Gln=3; end    %Khlebnikov 2018
        P.pool_names{1+P.n_cest_pools}='Gln_1';
        P.pool_names{2+P.n_cest_pools}='Gln_2';
        P.pool_names{3+P.n_cest_pools}='Gln_3';
        P.H=[P.H                1*P.Gln     1*P.Gln     3*P.Gln];    %3 exchangeable protons at pH=7
        P.kex=[P.kex            17          49          22880];      %Hz Khlebnikov 2018 at pH=7
        P.dw=[P.dw              2.15        2.87        3.18];         %ppm
        P.T1=[P.T1      P.T1_water  P.T1_water  P.T1_water]; %by default T1 is the same as water
        P.T2=[P.T2      0.014   0.014   0.014]; %Khlebnikov 2018
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=1*P.Gln;'));
        P=add_constrained_relation(P,strcat('P.H(',num2str(2+P.n_cest_pools),')=1*P.Gln;'));
        P=add_constrained_relation(P,strcat('P.H(',num2str(3+P.n_cest_pools),')=3*P.Gln;'));
        P.n_cest_pools=P.n_cest_pools+3; 
        
    case 'glucose_simplified_to_1pool' %Jin 2014
        if ~isfield(P,'Glc'); P.Glc=10; end    %if not given, init [Glc] at 10mM
        P.pool_names{1+P.n_cest_pools}='Glc_1';
        P.H=[P.H                5*P.Glc];    
        P.kex=[P.kex           6000];      %Hz   
        P.dw=[P.dw              1.5];         %ppm
        P.T1=[P.T1  P.T1_water  ]; 
        P.T2=[P.T2  0.008       ] ; 
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=5*P.Glc;'));
        P.n_cest_pools=P.n_cest_pools+1;
        
    case 'glucose_invitro'  %as given in Zaiss 2019, at approx pH=7 and 25° (phantoms)
        if ~isfield(P,'Glc'); P.Glc=10; end    %if not given, init [Glc] at 10mM
        if ~isfield(P,'AR'); P.AR=0.36; end    %if not given, init anomeric ratio at 0.36
        P.pool_names{1+P.n_cest_pools}='Glc_1';
        P.pool_names{2+P.n_cest_pools}='Glc_2';
        P.pool_names{3+P.n_cest_pools}='Glc_3';
        P.pool_names{4+P.n_cest_pools}='Glc_4';
        P.H=[P.H                P.Glc       3*P.Glc       P.AR*P.Glc    (1-P.AR)*P.Glc];    
        P.kex=[P.kex            2300        5000          3800          10000];      %Hz   
        P.dw=[P.dw              0.75        1.29          2.18          2.88];         %ppm %modified offset of 1st group compared to Zaiss (this fits results better)
        P.T1=[P.T1  P.T1_water  P.T1_water    P.T1_water    P.T1_water]; 
        P.T2=[P.T2  0.008       0.015         0.015         0.015] ; 
        
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=P.Glc;'));
        P=add_constrained_relation(P,strcat('P.H(',num2str(2+P.n_cest_pools),')=3*P.Glc;'));
        P=add_constrained_relation(P,strcat('P.H(',num2str(3+P.n_cest_pools),')=P.AR*P.Glc;'));
        P=add_constrained_relation(P,strcat('P.H(',num2str(4+P.n_cest_pools),')=(1-P.AR)*P.Glc;'));
        P.n_cest_pools=P.n_cest_pools+4;
        
    case 'glucose_physio'  %as given in Zaiss 2019 Table 2 (pH 7.2 and 37°C)
        if ~isfield(P,'Glc'); P.Glc=3; end    %if not given, init [Glc] at 3mM
        if ~isfield(P,'AR'); P.AR=0.38; end    %if not given, init anomeric ratio at 0.36
        P.pool_names{1+P.n_cest_pools}='Glc_1';
        P.pool_names{2+P.n_cest_pools}='Glc_2';
        P.pool_names{3+P.n_cest_pools}='Glc_3';
        P.pool_names{4+P.n_cest_pools}='Glc_4';
        P.H=[P.H                P.Glc       3*P.Glc       P.AR*P.Glc    (1-P.AR)*P.Glc];    
        P.kex=[P.kex            2900        6500          5200          14000];      %Hz   
        P.dw=[P.dw              0.74        1.29          2.18          2.88];         %ppm
        P.T1=[P.T1  P.T1_water  P.T1_water    P.T1_water    P.T1_water]; 
        P.T2=[P.T2  0.008       0.015         0.015         0.015] ; 
        
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=P.Glc;'));
        P=add_constrained_relation(P,strcat('P.H(',num2str(2+P.n_cest_pools),')=3*P.Glc;'));
        P=add_constrained_relation(P,strcat('P.H(',num2str(3+P.n_cest_pools),')=P.AR*P.Glc;'));
        P=add_constrained_relation(P,strcat('P.H(',num2str(4+P.n_cest_pools),')=(1-P.AR)*P.Glc;'));
        P.n_cest_pools=P.n_cest_pools+4;    

    case 'glucose_Khlebnikov'  %as given in Khlebnikov 2019 (PBS, pH=7, 37°C)
        if ~isfield(P,'Glc'); P.Glc=1.0; end    %if not given, init [Glc] at 1mM
        if ~isfield(P,'AR'); P.AR=0.43; end    %if not given, init anomeric ratio at 0.43
        P.pool_names{1+P.n_cest_pools}='Glc_1';
        P.pool_names{2+P.n_cest_pools}='Glc_2';
        P.pool_names{3+P.n_cest_pools}='Glc_3';
        P.pool_names{4+P.n_cest_pools}='Glc_4';
        P.pool_names{5+P.n_cest_pools}='Glc_5';
        P.H=[P.H                P.AR*P.Glc       (1-P.AR)*P.Glc       P.Glc    2*P.Glc  P.Glc];    
        P.kex=[P.kex            3860        4750          3940          2560    950];      %Hz   
        P.dw=[P.dw              2.18        2.88          1.10          1.39    0.74];         %ppm
        P.T1=[P.T1  P.T1_water  P.T1_water    P.T1_water    P.T1_water]; 
        P.T2=[P.T2  0.0069       0.0069         0.0069         0.0069   0.0069] ; 
        
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=P.AR*P.Glc;'));
        P=add_constrained_relation(P,strcat('P.H(',num2str(2+P.n_cest_pools),')=(1-P.AR)*P.Glc;'));
        P=add_constrained_relation(P,strcat('P.H(',num2str(3+P.n_cest_pools),')=P.Glc;'));
        P=add_constrained_relation(P,strcat('P.H(',num2str(4+P.n_cest_pools),')=2*P.Glc;'));
        P=add_constrained_relation(P,strcat('P.H(',num2str(5+P.n_cest_pools),')=P.Glc;'));
        P.n_cest_pools=P.n_cest_pools+5;   
        
    case 'glutamate_physio' %as given in Wermter 2015
        if ~isfield(P,'Glu'); P.Glu=10; end    %if not given, init [Glu] at 10mM    
        if ~isfield(P,'kex_Glu'); P.kex_Glu=7738; end
        if ~isfield(P,'T2_Glu'); P.T2_Glu=10; end
        
        P.pool_names{1+P.n_cest_pools}='Glu';
        P.H=[P.H                3*P.Glu];    
        P.kex=[P.kex            7738];      %Hz   
        P.dw=[P.dw              3.0];         %ppm
        P.T1=[P.T1  P.T1_water  ]; 
        P.T2=[P.T2  0.015       ] ; 
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=3*P.Glu;'));
        P=add_constrained_relation(P,strcat('P.T2(',num2str(1+P.n_cest_pools),')=P.T2_Glu;'));
        P=add_constrained_relation(P,strcat('P.kex(',num2str(1+P.n_cest_pools),')=P.kex_Glu;'));

        P.n_cest_pools=P.n_cest_pools+1;    

   case 'glutamate_physio_2' %as given in Wermter 2015
        if ~isfield(P,'Glu2'); P.Glu2=10; end    %if not given, init [Glu] at 10mM      
        
        P.pool_names{1+P.n_cest_pools}='Glu2';
        P.H=[P.H                3*P.Glu2];    
        P.kex=[P.kex            7738];      %Hz   
        P.dw=[P.dw              3.0];         %ppm
        P.T1=[P.T1  P.T1_water  ]; 
        P.T2=[P.T2  0.015       ] ; 
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=3*P.Glu2;'));
        
        P.n_cest_pools=P.n_cest_pools+1;   
        
    case 'glutamate_singlepool_pHdependent'
        if ~isfield(P,'Glu'); P.Glu=10; end    %if not given, init [Glu] at 10mM
        if ~isfield(P,'pH'); P.pH=7.0; end
        
        if ~isfield(P,'kex0_Glu'); P.kex0_Glu=3751; end
        if ~isfield(P,'kex1_Glu'); P.kex1_Glu=3987; end
        
        
        P.pool_names{1+P.n_cest_pools}='Glu';
        P.H=[P.H                3*P.Glu];    
        P.kex=[P.kex            P.kex0_Glu+P.kex1_Glu*10^(P.pH-7)];      %Hz   
        P.dw=[P.dw              3.0];         %ppm
        P.T1=[P.T1  P.T1_water  ]; 
        P.T2=[P.T2  0.015       ] ; 
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=3*P.Glu;'));
        P=add_constrained_relation(P,strcat('P.kex(',num2str(1+P.n_cest_pools),')=P.kex0_Glu+P.kex1_Glu*10^(P.pH-7);'));
        
        P.n_cest_pools=P.n_cest_pools+1;
    
    case 'glutamate_OH_and_NH3_pHdependent'
        if ~isfield(P,'Glu'); P.Glu=10; end    %if not given, init [Glu] at 10mM
        if ~isfield(P,'pH'); P.pH=7.0; end
        if ~isfield(P,'pKa2'); P.pKa2=4.25; end
        
%         if ~isfield(P,'kex0_Glu_1'); P.kex0_Glu_1=1800; end
%         if ~isfield(P,'kex1_Glu_1'); P.kex1_Glu_1=3000; end

        if ~isfield(P,'kex0_Glu'); P.kex0_Glu_1=2100; end
        if ~isfield(P,'kex1_Glu'); P.kex1_Glu_1=3900; end
        
        
        P.pool_names{1+P.n_cest_pools}='Glu';
        P.pool_names{2+P.n_cest_pools}='Glu_OH';
        P.H=[P.H                3*P.Glu   P.Glu*10^(P.pKa2-P.pH)/(1+10^(P.pKa2-P.pH))];    
        P.kex=[P.kex            P.kex0_Glu_1+P.kex1_Glu_1*10^(P.pH-7)           2500];      %Hz   
        P.dw=[P.dw              3.0             1.0];         %ppm
        P.T1=[P.T1  P.T1_water  P.T1_water   ]; 
        P.T2=[P.T2  0.015         0.015  ] ; 
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=3*P.Glu;'));
        P=add_constrained_relation(P,strcat('P.kex(',num2str(1+P.n_cest_pools),')=P.kex0_Glu+P.kex1_Glu*10^(P.pH-7);'));
        P=add_constrained_relation(P,strcat('P.H(',num2str(2+P.n_cest_pools),')=P.Glu*10^(P.pKa2-P.pH)/(1+10^(P.pKa2-P.pH));'));
        
        P.n_cest_pools=P.n_cest_pools+2;

    case 'vesicularGlu'
        if ~isfield(P,'fves'); P.fves=0.001; end
        if ~isfield(P,'kex_ves'); P.kex_ves=300; end
        if ~isfield(P,'T2_ves'); P.T2_ves=0.015; end

        P.pool_names{1+P.n_cest_pools}='Glu_ves';
        P.H=[P.H                P.fves*P.H_water  ];    
        P.kex=[P.kex            P.kex_ves];      %Hz   
        P.dw=[P.dw              3.0      ];         %ppm
        P.T1=[P.T1  P.T1_water  ]; 
        P.T2=[P.T2  0.005       ] ; 
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=P.fves*P.H_water;'));
        P=add_constrained_relation(P,strcat('P.kex(',num2str(1+P.n_cest_pools),')=P.kex_ves;'));
        P=add_constrained_relation(P,strcat('P.T2(',num2str(1+P.n_cest_pools),')=P.T2_ves;'));
        
        P.n_cest_pools=P.n_cest_pools+1;

   case 'LipoGlu'
        if ~isfield(P,'fves'); P.fves=0.001; end
        if ~isfield(P,'kex_ves'); P.kex_ves=300; end
        if ~isfield(P,'T2_ves'); P.T2_ves=0.01; end

        if ~isfield(P,'Glu_ves'); P.Glu_inves=100; end
        if ~isfield(P,'kex_Glu_ves'); P.kex_Glu_ves=1800; end
        if ~isfield(P,'T2_Glu_ves'); P.T2_Glu_ves=0.01; end

        P.intramol_transfer_modeling=0;
        P.analytic=0; P.numeric=1;
        P.index_lipoGlu=1+P.n_cest_pools;
        % P.kex_intramol=[0, P.kex_Glu_ves;
        %     (3*P.Glu_inves/(2*55556))*P.kex_Glu_ves, 0];
                        

        P.pool_names{1+P.n_cest_pools}='H2O_ves';
        P.pool_names{2+P.n_cest_pools}='Glu_inves';
        
        P.H=[P.H                P.fves*P.H_water  3*P.Glu_inves];    
        P.kex=[P.kex            P.kex_ves   0];      %Hz   
        P.dw=[P.dw              0.0     3.0 ];         %ppm
        P.T1=[P.T1  P.T1_water  P.T1_water]; 
        P.T2=[P.T2  P.T2_ves       P.T2_Glu_ves] ; 
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=P.fves*P.H_water;'));
        P=add_constrained_relation(P,strcat('P.kex(',num2str(1+P.n_cest_pools),')=P.kex_ves;'));
        P=add_constrained_relation(P,strcat('P.T2(',num2str(1+P.n_cest_pools),')=P.T2_ves;'));

        P=add_constrained_relation(P,strcat('P.H(',num2str(2+P.n_cest_pools),')=3*P.Glu_inves;'));
        P=add_constrained_relation(P,strcat('P.kex(',num2str(2+P.n_cest_pools),')=0;'));
        P=add_constrained_relation(P,strcat('P.T2(',num2str(2+P.n_cest_pools),')=P.T2_Glu_ves;'));

        %intramolecular exchanges
        P=add_constrained_relation(P,'P.kex_intramol=[0, P.kex_Glu_ves;(3*P.Glu_inves/(2*55556))*P.kex_Glu_ves, 0];');
        
        %constrain vesicular water T1/T2 to the same values as bulk water?
        P=add_constrained_relation(P,strcat('P.T2_ves=P.T2_water;'));
        P=add_constrained_relation(P,strcat('P.T1(',num2str(1+P.n_cest_pools),')=P.T1_water;'));

        P.n_cest_pools=P.n_cest_pools+2;

    case 'UnknownPool'
        if ~isfield(P,'H_Unkn'); P.H_Unkn=5; end          
        
        P.pool_names{1+P.n_cest_pools}='Unkn';
        P.H=[P.H                P.H_Unkn];    
        P.kex=[P.kex            5000];      %Hz   
        P.dw=[P.dw              3.0];         %ppm
        P.T1=[P.T1  P.T1_water  ]; 
        P.T2=[P.T2  0.010       ] ; 
        
        P=add_constrained_relation(P,strcat('P.H(',num2str(1+P.n_cest_pools),')=P.H_Unkn;'));
        
        P.n_cest_pools=P.n_cest_pools+1;   

    otherwise
        fprintf('Unknown metabolite model');
end