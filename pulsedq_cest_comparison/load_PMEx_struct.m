function PMEX=load_PMEx_struct(P,seq_name,sim_name)

seqid = 'Test_Cecile';
seq_defs.n_pulses      = P.n              ; % number of pulses
seq_defs.tp            = P.tp           ; % pulse duration [s]
seq_defs.td            = P.td            ; % interpulse delay [s]
seq_defs.Trec          = P.Trec             ; % recovery time [s]
seq_defs.Trec_M0       = P.Trec             ; % recovery time before M0 [s]
seq_defs.M0_offset     = P.normalized          ; % m0 offset [ppm]
seq_defs.DCsat         = (seq_defs.tp)/(seq_defs.tp+seq_defs.td); % duty cycle
seq_defs.offsets_ppm   = P.xZspec; % offset vector [ppm]
seq_defs.num_meas      = numel(seq_defs.offsets_ppm)   ; % number of repetition
seq_defs.Tsat          = seq_defs.n_pulses*(seq_defs.tp+seq_defs.td) - ...
                         seq_defs.td ;  % saturation time [s]
seq_defs.B0            = P.FREQ/gamma_               ; % B0 [T]
seq_defs.seq_id_string = seqid           ; % unique seq id

lims = getScannerLimits();

B1pa          = P.B1;                  % mean sat pulse b1 [uT]
gyroRatio_hz  = 42.5764;               % for H [Hz/uT]
gyroRatio_rad = gyroRatio_hz*2*pi;     % [rad/uT]
fa_sat        = B1pa*gyroRatio_rad*seq_defs.tp; % flip angle of sat pulse
% satPulse      = mr.makeSincPulse(fa_sat, 'Duration', seq_defs.tp, 'system', lims,'timeBwProduct', 2,'apodization', 0.15);
satPulse      = mr.makeBlockPulse(fa_sat, 'Duration', seq_defs.tp, 'system', lims);

[B1cwpe,B1cwae,B1cwae_pure,alpha]= calculatePowerEquivalents(satPulse,seq_defs.tp,seq_defs.td,0,gyroRatio_hz);

seq_defs.B1cwpe = B1cwpe;

% init sequence with system limits
seq = SequenceSBB(lims);

% pulseq uses offsets in Hz
offsets_Hz = seq_defs.offsets_ppm*gyroRatio_hz*seq_defs.B0;

% loop through offsets
for currentOffset = offsets_Hz
    % if the current offset is an M0 scan we use the M0 recovery time
    if currentOffset == seq_defs.M0_offset*gyroRatio_hz*seq_defs.B0
        if seq_defs.Trec_M0 > 0
            seq.addBlock(mr.makeDelay(seq_defs.Trec_M0));
        end
    else % if no M0 scan use Trec
        if seq_defs.Trec > 0
            seq.addBlock(mr.makeDelay(seq_defs.Trec)); % recovery time
        end
    end
    % set frequency offset of the pulse
    satPulse.freqOffset = currentOffset;
    
    % take care of the accumulated phase during the saturation
    accumPhase=0;
    
    % loop through pulses
    for np = 1:seq_defs.n_pulses
        % set accumulated pahse from previous rf pulse
        satPulse.phaseOffset = mod(accumPhase,2*pi); 
        % add satURATION PULSE pulse
        seq.addBlock(satPulse) 
        % calculate the phase for the next rf pulse
        accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6),2*pi);
        % delay between pulses
        if np < seq_defs.n_pulses
            seq.addBlock(mr.makeDelay(seq_defs.td)); % add delay
        end
    end
    % add the standard spoiling gradients
    seq.addSpoilerGradients();
    % add the readout trigger event
    seq.addPseudoADCBlock(); 
end

%% set definitions
def_fields = fieldnames(seq_defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, seq_defs.(def_fields{n_id}));
end
% write the pulseq-file with author info
seq_filename = 'Test_sequence.seq';
author = 'CecileMaguin';
seq.write(seq_filename, author);

param_filename = 'C:/Users/cm268688/Desktop/KHerz_pulsedqCEST/pulseq-cest-library/sim-library/Test_Cecile_glc.yaml';
% sim_parameters = yaml.ReadYaml(param_filename);
% 
% 
% %water pool 
% sim_parameters.water_pool.t1=P.T1_water;
% sim_parameters.water_pool.t2=P.T2_water;
% 
% for i=1:P.n_cest_pools
%     sim_parameters.cest_pool.(P.pool_names{i}).f=P.CALC.fH(i);
%     sim_parameters.cest_pool.(P.pool_names{i}).t1=P.T1(i);
%     sim_parameters.cest_pool.(P.pool_names{i}).t2=P.T2(i);
%     sim_parameters.cest_pool.(P.pool_names{i}).k=P.kex(i);
%     sim_parameters.cest_pool.(P.pool_names{i}).k=P.dw(i);
% end

PMEX = readSimulationParameters(param_fn);

end

