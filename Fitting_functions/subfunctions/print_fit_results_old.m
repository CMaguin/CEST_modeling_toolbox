function print_fit_results_old(Results,OPTIM, vars)

    if strcmp(vars,'all')
        
        %Global vars
        for i=1:numel(Results.global_vars)
            fprintf(strcat(Results.global_vars{i},'= %.2f +- %.2f ', get_unit(Results.global_vars{i}),' \n'),Results.(Results.global_vars{i}));
        end
        
        %Voxel vars mean values
        for i=1:numel(Results.global_vars)
            fprintf('[Glu]= %.2f mM +- %.2f mM (Original= %.2f mM) \n',Results.voxel_val(:,strcmp(Results.voxel_vars,'Glu')), Results.ci_voxel(:,strcmp(Results.voxel_vars,'Glu'),2), str2double(cell2mat(OPTIM.voxel_start(strcmp(OPTIM.voxel_vars,'Glu')))) );
        end
    else
        try
            fprintf('Avg %s = %.2f %s +- %.2f mM (Original= %.2f mM) \n',vars, mean(Results.global_val(:,strcmp(Results.global_vars,vars))), get_unit(vars), mean(Results.ci_global(:,strcmp(Results.global_vars,vars),2)), str2double(cell2mat(OPTIM.global_start(strcmp(OPTIM.global_vars,vars)))) );
        catch
            try
                fprintf('Avg %s = %.2f %s +- %.2f mM (Original= %.2f mM) \n',vars, mean(Results.voxel_val(:,strcmp(Results.voxel_vars,vars))), get_unit(vars), mean(Results.ci_voxel(:,strcmp(Results.voxel_vars,vars),2)), str2double(cell2mat(OPTIM.voxel_start(strcmp(OPTIM.voxel_vars,vars)))) );
            catch
                fprintf('Avg %s = %.2f %s +- %.2f mM (Original= %.2f mM) \n',vars, mean(Results.indiv_val(:,strcmp(Results.indiv_vars,vars))), get_unit(vars), mean(Results.ci_indiv(:,strcmp(Results.indiv_vars,vars),2)), str2double(cell2mat(OPTIM.indiv_start(strcmp(OPTIM.indiv_vars,vars)))) );
            end
        end

end