function print_fit_results(Results,OPTIM,P, vars)

    if strcmp(vars,'all')
        
        %Global vars
        for i=1:numel(OPTIM.global_vars)
            fprintf(strcat(OPTIM.global_vars{i},'= %.2f +- %.2f ', get_unit(OPTIM.global_vars{i}),' \n'),Results.(OPTIM.global_vars{i}));
        end
        
        %Voxel vars mean values
        for i=1:numel(Results.voxel_vars)
            fprintf('Mean %s= %.2f  +- %.2f %s (Original= %.2f  %s , LB= %.2f %s, UB=%.2f %s) \n',Results.voxel_vars{i},mean(Results.voxel_val(:,strcmp(Results.voxel_vars,Results.voxel_vars{i})),'all'), mean(Results.ci_voxel(:,strcmp(Results.voxel_vars,Results.voxel_vars{i}),2),'all'), get_unit(Results.voxel_vars{i}),str2double(cell2mat(OPTIM.voxel_start(strcmp(OPTIM.voxel_vars,Results.voxel_vars{i})))), get_unit(Results.voxel_vars{i}),str2double(cell2mat(OPTIM.voxel_lb(strcmp(OPTIM.voxel_vars,Results.voxel_vars{i})))), get_unit(Results.voxel_vars{i}),str2double(cell2mat(OPTIM.voxel_ub(strcmp(OPTIM.voxel_vars,Results.voxel_vars{i})))), get_unit(Results.voxel_vars{i}) );
        end
    else
        if isfield(Results, vars)
            fprintf('Avg %s = %.2f +- %.2f %s (Original= %.2f %s) \n',vars, Results.vars,  Results.ci(:,strcmp(OPTIM.global_vars,vars),2), get_unit(vars),str2double(cell2mat(OPTIM.global_start(strcmp(OPTIM.global_vars,vars)))) );
        else
            if any(strcmp(Results.voxel_vars,vars))
                fprintf('Avg %s = %.2f  +- %.2f %s (Original= %.2f ) \n',vars, mean(Results.voxel_val(:,strcmp(Results.voxel_vars,vars)),'all'),  mean(Results.ci_voxel(:,strcmp(Results.voxel_vars,vars),2),'all'), get_unit(vars),str2double(cell2mat(OPTIM.voxel_start(strcmp(OPTIM.voxel_vars,vars)))) );
            else
                if any(strcmp(Results.indiv_vars,vars))
                    fprintf('Avg %s = %.2f  +- %.2f %s (Original= %.2f ) \n',vars, mean(Results.indiv_val(:,strcmp(Results.indiv_vars,vars)),'all'), mean(Results.ci_indiv(:,strcmp(Results.indiv_vars,vars),2),'all'),get_unit(vars), str2double(cell2mat(OPTIM.indiv_start(strcmp(OPTIM.indiv_vars,vars)))) );
                else
                    fprintf('%s was fixed at %.2f %s \n', vars, P.(vars),get_unit(vars))
                end
            end
        end

    end