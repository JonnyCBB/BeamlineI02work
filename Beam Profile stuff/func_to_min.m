%% Function to calculate the rmsd of calculated values from a gaussain with experimental values using given parameters
function val_to_min = func_to_min(noise,data)
%experimental data
beam_pos = data(:,1);
flux_data = data(:,2);

%% Create convolution variable - corresponds to aperature diameter and translation

%Aperature diameter
ap_diam = 0.01;

%Aperature translation
ap_trans = beam_pos(2) - beam_pos(1);

%Coverage of aperature
cov_of_ap = ap_diam/ap_trans;
%Convert the value from a double class to an integer
cov_of_ap = uint8(cov_of_ap);

%Convolution variable - i.e the aperature
aperature = ones(cov_of_ap,1);

% aperature = ones(3,1);

%% Deconvolve the data signal with the aperature
sig = deconvwnr(flux_data,aperature,noise);

%% Calculate the Gaussian parameters for the best fit curve

% %Only use if using the mygaussfit function
% [coeff_wie, mu_wie, sigma_wie] = mygaussfit(beam_pos, sig);

%Only use if using the gaussfit function
[sigma_wie, mu_wie] = gaussfit(beam_pos, sig);

%% Calculate the Gaussian function (making sure we scale the Gaussian correctly)

% %Only use if using the mygaussfit function
% flux_wie = coeff_wie * exp( - (beam_pos-mu_wie).^2 / (2*sigma_wie^2));

%Only use if using the gaussfit function
%Gaussian curve
flux_wie = 1/(sqrt(2*pi)* sigma_wie ) * exp( - (beam_pos-mu_wie).^2 / (2*sigma_wie^2));
%Compute scale factor
scale_fac = max(sig)/max(flux_wie);
%Scale the gaussian
flux_wie = scale_fac*flux_wie;

%% Convolve the calculated flux signal with the aperature
calc_obs_sig = conv(flux_wie,aperature,'same');

%% Find RMSD
val_to_min = rmsd(calc_obs_sig,flux_data);

end

