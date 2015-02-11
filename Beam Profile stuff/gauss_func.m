%% Gaussian Function
function gauss_val = gauss_func(x,param)
gauss_val = param(1).*exp(-((x-param(2)).^2)./(2*(param(3).^2)));
end

