%% Script to Test which distribution fits the x-ray beam profile the best
clear all
close all
clc

%% User Input

%Take initial guess for the parameters for the 
height_of_peak = 0.035;
x_location_of_peak = 3.56;
FWHM = 0.07;

height_of_peak = 0.015;
x_location_of_peak = 2.55;
FWHM = 0.07;

height_of_peak = 0.06;
x_location_of_peak = 50.05;
FWHM = 0.01;
% % 
% % height_of_peak = 0.035;
% % x_location_of_peak = 1.1;
% % FWHM = 0.02;

%Calculate the scale parameters from the FWHM
gamma_lorentzian = FWHM/2;
standard_deviation_gaussian = FWHM/(2*sqrt(2*log(2)));

%form the initial guesses for the functions
x0_lorentzian = [height_of_peak; x_location_of_peak; gamma_lorentzian];
x0_gaussian = [height_of_peak; x_location_of_peak; standard_deviation_gaussian ];
%% Read in the beam profile data

%Read x data file to import the data 
importfile('13765.dat');

%Create relevant variables from the data file
%make variable for the postion of the beam
beam_pos = data(:,1);

%make variable for the flux data
flux_data = data(:,2);
if min(flux_data) < 0 
    num_to_add_to_data = abs(min(flux_data));
    flux_data = flux_data + num_to_add_to_data;
end

%% Calculate the least squares fits

%Calculate the least squares fit for the Lorentzian
[x_lorentzian,resnorm_lorentzian] = lsqcurvefit(@lorentzian_func,x0_lorentzian,beam_pos,flux_data);
%Calculate the least squares fit for the Gaussian
[x_gaussian,resnorm_gaussian] = lsqcurvefit(@gaussian_func,x0_gaussian,beam_pos,flux_data);

%% Plot the functions
F_lorentz = lorentzian_func(x_lorentzian,beam_pos);
F_gauss = gaussian_func(x_gaussian,beam_pos);

%Plot the data for the y direction
figure('name','Distribution Fits to data')
plot(beam_pos,flux_data,'o',beam_pos,F_lorentz,'r',beam_pos,F_gauss,'g')
title('Figure showing the fits of the distributions to the data','FontSize',14)
ylabel('ipin reading (Arbitrary Units)','FontSize',14)
xlabel('Position (Millimetres)','FontSize',14)
h2 = legend('miniap reading','Fitted Lorentzian','Fitted Gaussian');
set(h2,'FontSize',12)

%% Print to screen
FWHM_G = 2*x_gaussian(3)*sqrt(2*log(2));
FWHM_L = 2*x_lorentzian(3);

fprintf('The FWHM according to the Gaussian is %f\n',FWHM_G)
fprintf('The FWHM according to the Lorentzian is %f\n\n',FWHM_L)

fprintf('The squared 2-norm of residuals for Lorentzian is %f\n',resnorm_lorentzian)
fprintf('The squared 2-norm of residuals for Gaussian is %f\n\n',resnorm_gaussian)

if resnorm_lorentzian < resnorm_gaussian
    fprintf('The Lorentzian is a better fit than the Gaussian\n\n')
elseif resnorm_lorentzian > resnorm_gaussian
    fprintf('The Gaussian is a better fit than the Lorentzian\n\n')
else
    fprintf('Both fits are equally as good\n\n')
end

FWHM_G = 2*x_gaussian(3)*sqrt(2*log(2));
FWHM_L = 2*x_lorentzian(3);