function [flux_inv beam_pos] = beam_FWHM_Wiener(filename)
%% Parameter to optimise: Take an initial guess
noise_param_guess = 0.7;

%% Read select file to import the data 
full_data = importdata(filename);

%% Create relevant variables from the data file
%make variable for the postion of the beam
beam_pos = full_data.data(:,1);

%make variable for the flux data
flux_data = full_data.data(:,2);

%% Calculate the Optimal Noise Parameter
%This function minimises the Root Mean Squared Deviation of the original
%flux data from the .dat file with the values that are calculated from
%ideal gaussian distributions.
[noise_param rmsd]= fminsearch(@(x) func_to_min(x,full_data.data),noise_param_guess);

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

%% Deconvolve the data signal with the aperature

[signal rem] = deconv(flux_data,aperature);
sig = deconvwnr(flux_data,aperature,noise_param);

%% Calculate the Gaussian parameters for the best fit curve

[sigma_wie, mu_wie] = gaussfit(beam_pos, sig);

[sigma_inv, mu_inv] = gaussfit(beam_pos(3:end-2), signal);

%% Calculate the Gaussian function (making sure we scale the Gaussian correctly)

%Gaussian curve
flux_wie = 1/(sqrt(2*pi)* sigma_wie ) * exp( - (beam_pos-mu_wie).^2 / (2*sigma_wie^2));
%Compute scale factor
scale_fac = max(sig)/max(flux_wie);
%Scale the gaussian
flux_wie = scale_fac*flux_wie;

%Calculate the Gaussian Curve
flux_inv = 1/(sqrt(2*pi)* sigma_inv ) * exp( - (beam_pos-mu_inv).^2 / (2*sigma_inv^2));
%Compute Scale factor
scale_fac = max(signal)/max(flux_inv);
%Scale the Gaussian curve
flux_inv = scale_fac*flux_inv;

%% Convolute the original signal with the aperature

calc_obs_signal = conv(flux_inv,aperature,'same');
calc_obs_sig = conv(flux_wie,aperature,'same');

%% Fit guassian to the colvolved signals

[sigma_wie_con, mu_wie_con] = gaussfit(beam_pos, calc_obs_sig);
[sigma_inv_con, mu_inv_con] = gaussfit(beam_pos, calc_obs_signal);

%% Calculate the Gaussian function for convolved signal (making sure we scale the Gaussian correctly)

%Gaussian curve
con_flux_wie = 1/(sqrt(2*pi)* sigma_wie_con ) * exp( - (beam_pos-mu_wie_con).^2 / (2*sigma_wie_con^2));
%Compute scale factor
scale_fac_con = max(calc_obs_sig)/max(con_flux_wie);
%Scale the gaussian
con_flux_wie = scale_fac_con*con_flux_wie;

%Calculate the Gaussian Curve
con_flux_inv = 1/(sqrt(2*pi)* sigma_inv_con ) * exp( - (beam_pos-mu_inv_con).^2 / (2*sigma_inv_con^2));
%Compute Scale factor
scale_fac_con = max(calc_obs_signal)/max(con_flux_inv);
%Scale the Gaussian curve
con_flux_inv = scale_fac_con*con_flux_inv;

%% Calculate the FWHM for the Gaussians
% fprintf( 'Data is not normalized! The pdf sums to: %f. Normalizing...\n\r', s );

FWHM_flux_inv = 2*sigma_inv*sqrt(2*log(2))
FWHM_flux_inv_con = 2*sigma_inv_con*sqrt(2*log(2))

FWHM_flux_wie = 2*sigma_wie*sqrt(2*log(2))
FWHM_flux_wie_con = 2*sigma_wie_con*sqrt(2*log(2))

%% Plot All of the Data in the Figures

figure('name','Signal Plots for wiener deconvolution')
plot(beam_pos,flux_data,'o',beam_pos,sig,'x',beam_pos,flux_wie,'r',beam_pos,calc_obs_sig,'k')
title('Amp reading against aperature Position','FontSize',14)
ylabel('ipin reading (microamps)','FontSize',14)
xlabel('Aperature position (millimetres)','FontSize',14)
h = legend('original data','deconvoluted signal','gaussian fit (deconv sig)','Convolution gaussian fit with aperature');
set(h,'FontSize',12)

figure('name','Signal Plots')
plot(beam_pos,flux_data,'o',beam_pos(3:end-2),signal,'x',beam_pos,flux_inv,'r',beam_pos,calc_obs_signal,'k')
title('Amp reading against aperature Position','FontSize',14)
ylabel('ipin reading (microamps)','FontSize',14)
xlabel('Aperature position (millimetres)','FontSize',14)
h = legend('original data','deconvoluted signal','gaussian fit (deconv sig)','Convolution gaussian fit with aperature');
set(h,'FontSize',12)
%% Calculate the FWHM of fit with the original data

[sigma, mu] = gaussfit(beam_pos, flux_data);
FWHM = 2*sigma*sqrt(2*log(2))

%% Calculate the cumulative sum of the deconvoluted data

cum_flux_wie = cumsum(flux_wie);
cum_flux_inv = cumsum(flux_inv);

%% Plot the cumulative sum flux

figure('name','Cumulative Flux Plots')
plot(beam_pos,cum_flux_wie,'b',beam_pos,cum_flux_inv,'r')
title('Cumulative Flux plots against the beam position')
ylabel('Cumulative i_pin reading')
xlabel('Beam Position')
legend('wiener deconv flux','inverse deconv flux')
end

