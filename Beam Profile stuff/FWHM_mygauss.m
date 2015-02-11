%% Script for deconvolution of beam data to obtain true FWHM
close all
clear all
clc

%% Parameter to optimise: Take an initial guess
noise_param_guess = 0.7;

%% Read select file to import the data 
importfile('beam_prof_x.dat');

%% Create relevant variables from the data file
%make variable for the postion of the beam
beam_pos = data(:,1);

%make variable for the flux data
flux_data = data(:,2);

%% Calculate the Optimal Noise Parameter
%This function minimises the Root Mean Squared Deviation of the original
%flux data from the .dat file with the values that are calculated from
%ideal gaussian distributions.
[noise_param rmsd]= fminsearch(@(x) func_to_min(x,data),noise_param_guess);

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

[coeff_wie, mu_wie, sigma_wie] = mygaussfit(beam_pos, sig);

[coeff_inv, mu_inv, sigma_inv] = mygaussfit(beam_pos(3:end-2), signal);

%% Calculate the Gaussian function (making sure we scale the Gaussian correctly)

%Gaussian curve
flux_wie = coeff_wie * exp( - (beam_pos-mu_wie).^2 / (2*sigma_wie^2));

%Calculate the Gaussian Curve
flux_inv = coeff_inv * exp( - (beam_pos-mu_inv).^2 / (2*sigma_inv^2));

% %% Plot the beam profiles
% 
% figure('name','Signal Plots for wiener deconvolution')
% plot(beam_pos,flux_data,'o',beam_pos,sig,'x',beam_pos,flux_wie,'r')
% 
% figure('name','Signal Plots')
% plot(beam_pos,flux_data,'o',beam_pos(3:end-2),signal,'x',beam_pos,flux_inv,'r')

%% Convolute the original signal with the aperature

calc_obs_signal = conv(flux_inv,aperature,'same');
calc_obs_sig = conv(flux_wie,aperature,'same');

%% Fit guassian to the colvolved signals

[coeff_wie_con, mu_wie_con, sigma_wie_con] = mygaussfit(beam_pos, calc_obs_sig);

[coeff_inv_con, mu_inv_con, sigma_inv_con] = mygaussfit(beam_pos, calc_obs_signal);

%% Calculate the Gaussian function for convolved signal (making sure we scale the Gaussian correctly)

%Gaussian curve
con_flux_wie = coeff_wie_con * exp( - (beam_pos-mu_wie_con).^2 / (2*sigma_wie_con^2));

%Calculate the Gaussian Curve
con_flux_inv = coeff_inv_con * exp( - (beam_pos-mu_inv_con).^2 / (2*sigma_inv_con^2));

%% Calculate the FWHM for the Gaussians

FWHM_flux_inv = 2*sigma_inv*sqrt(2*log(2))
FWHM_flux_inv_con = 2*sigma_inv_con*sqrt(2*log(2))

FWHM_flux_wie = 2*sigma_wie*sqrt(2*log(2))
FWHM_flux_wie_con = 2*sigma_wie_con*sqrt(2*log(2))

%% Plot All of the Data in the Figures

figure('name','Signal Plots for wiener deconvolution')
plot(beam_pos,flux_data,'o',beam_pos,sig,'x',beam_pos,flux_wie,'r',beam_pos,calc_obs_sig,'k')
title('Amp reading against aperature Position')
ylabel('ipin reading (microamps)')
xlabel('Aperature position (millimetres)')
legend('original data','deconvoluted signal','gaussian fit (deconv sig)','Convolution gaussian fit with aperature')

figure('name','Signal Plots')
plot(beam_pos,flux_data,'o',beam_pos(3:end-2),signal,'x',beam_pos,flux_inv,'r',beam_pos,calc_obs_signal,'k')
title('Amp reading against aperature Position')
ylabel('ipin reading (microamps)')
xlabel('Aperature position (millimetres)')
legend('original data','deconvoluted signal','gaussian fit (deconv sig)','Convolution gaussian fit with aperature')

%% Calculate the FWHM of fit with the original data

[coeff, mu, sigma] = mygaussfit(beam_pos, flux_data);
FWHM = 2*sigma*sqrt(2*log(2))