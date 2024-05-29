function unit=get_unit(vars)
if numel(vars)>3
    if strcmp(vars(1:3),'kex')
        unit='Hz';
    elseif strcmp(vars(1:2),'dw')
        unit='ppm';
    elseif strcmp(vars(1:2),'H_')
        unit='mM';
    elseif strcmp(vars(1:3),'T1_') || strcmp(vars(1:3),'T2_')
        unit='s';
    else 
        unit='';
    end
elseif strcmp(vars,'Glc') || strcmp(vars,'Glu') || strcmp(vars,'Cr') || strcmp(vars,'PCr') || strcmp(vars,'APT') || strcmp(vars,'MI')
    unit='mM';
elseif strcmp(vars,'B1')
    unit='uT';
elseif strcmp(vars,'tp')    
    unit='s';
elseif strcmp(vars,'T1_water') || strcmp(vars,'T2_water')
    unit='s';
else
    unit='';
end
    

