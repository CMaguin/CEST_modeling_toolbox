function y=normc(x)
    if ~ismatrix(x)
        error('Input argument must be a 2D matrix');
    end
    y=zeros(size(x));
    scales = 1./sqrt(sum(abs(x).^2));
    for i=1:numel(scales)
        y(:,i)=x(:,i)*scales(i);
    end
end