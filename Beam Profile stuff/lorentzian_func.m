%% Function that calculates the 3 parameter Lorenztian function.
%x(1) = I (height of the peak)
%x(2) = x0 (location of peak)
%x(3) = gamma (scale parameter which is the half-width at half maximum (HWHM))
function F = lorentzian_func(x,xdata)
F = x(1).*(((x(3))^2)./((xdata - x(2)).^2 + (x(3))^2));
end

