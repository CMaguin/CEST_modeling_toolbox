function c=pc2c(parameter, pool_name,P)
if strcmp(pool_name,'Water')
    c=strcat(parameter,'_water');
elseif strcmp(pool_name,'MT')
    c=strcat(parameter,'_MT');
else
    i=find(strcmp(P.pool_names, pool_name));
    if isempty(i); warning('Unknown pool name'); end
    
    c=strcat(parameter, '(',num2str(i),')');

end

