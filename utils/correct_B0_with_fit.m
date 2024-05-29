function [Z, dw_fit]=correct_B0_with_fit(Z,xZspec, P, OPTIM)

if nargin<4
    if size(Z)<3
        OPTIM.nvoxels=1;
        [OPTIM.nconditions, ~]=size(Z);
    else
        [OPTIM.nvoxels, OPTIM.nconditions, ~]=size(Z);
    end
    OPTIM.vary_voxel={}; OPTIM.varyval_voxel=[];OPTIM.vary_indiv={};OPTIM.varyval_indiv=[];
end
if nargin<3
    P=init_Sim_struct();
    P.xZspec=xZspec;
    P.FREQ          = 11.7*gamma_;       % frequency (=B0[T] * gamma)
    P.B1            = 5;                  % B1 value in Ã‚µT
    P.TR=5;
    P.Trec          = 4;                    % recovery time in s
    P.Zi            = 1.0;                    % initial magnetisation (should be between -1 and +1)
    P.M0corr        = 1.0;                    % correction on M0

    P.pulsed  = 0;                    % 0 = cw saturation, 1 = pulsed saturation
    P.tp    = 1;                       % saturation time in s
    P.n     = 1;                        % choose n=1 for cw saturation
    P.shape = 'block';                  % choose 'block' for cw saturation (cases: SPINLOCK, seq_gauss, block, block_trap, gauss, sech, sinc_1, sinc_2, sinc_3, sinc_4)
    P.td=0 ;   
end

OPTIM_correction=init_OPTIM_struc(OPTIM.nvoxels, OPTIM.nconditions);
OPTIM_correction.fit_type=2;
OPTIM_correction.rescalefactMTR=10;

OPTIM_correction.vary_voxel=OPTIM.vary_voxel;
OPTIM_correction.varyval_voxel=OPTIM.varyval_voxel;
OPTIM_correction.vary_indiv=OPTIM.vary_indiv;
OPTIM_correction.varyval_indiv=OPTIM.varyval_indiv;

OPTIM_correction=add_voxel_var(OPTIM_correction, 'dw_water','0','-1.5','1.5');
OPTIM_correction=define_start_and_bounds(OPTIM_correction,P);


OPTIM_correction.fit_type=2; %take MTRasym to have better B0 inhomogeneity inferring
OPTIM_correction.rescalefactMTR=10; 

P.sim_type=1;

[Results,~, ~]=fit_data(P,Z,OPTIM_correction);
dw_fit=Results.voxel_val;

for voxel=1:OPTIM.nvoxels
   dwvox=Results.voxel_val(voxel);
   Z=correct_B0(Z,xZspec, dwvox);
end