%% Function that calculates the 3 parameter Gaussian function.
%x(1) = A (height of the peak)
%x(2) = mu (location of peak - The mean)
%x(3) = sigma (Standard Deviation)
function F = gaussian_func(x,xdata)
F = x(1).*exp(-((xdata - x(2)).^2)./(2*(x(3))^2));
end

