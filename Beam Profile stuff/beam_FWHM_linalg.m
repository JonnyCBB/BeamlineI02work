function [flux_real beam_pos, FWHM_real] = beam_FWHM_linalg(filename,gon_data)
%% Read select file to import the data 
full_data = importdata(filename);

%% Create relevant variables from the data file
%make variable for the postion of the beam
beam_pos = full_data.data(:,1);

%make variable for the flux data
flux_data = full_data.data(:,2);

%% Create aperture data

%Aperture diameter
ap_diam = 0.01;

%Aperture translation
ap_trans = beam_pos(2) - beam_pos(1);

%Coverage of aperture
cov_of_ap = ap_diam/ap_trans;
%Convert the value from a double class to an integer
cov_of_ap = uint8(cov_of_ap);

%% Fit Gaussian to the original flux data to obtain the Gaussian parameters

[sigma, mu] = gaussfit(beam_pos, flux_data);

%% Calculate the Gaussian function (making sure we scale the Gaussian correctly)

%Gaussian curve
flux = 1/(sqrt(2*pi)* sigma ) * exp( - (beam_pos-mu).^2 / (2*sigma^2));
%Compute scale factor
scale = max(flux_data)/max(flux);
%Scale the gaussian
flux = scale*flux;

% %Fit Gaussian using the least squares fit
% flux = amplitude * exp( - (beam_pos-mu).^2 / (2*sigma^2));

%% Set up matrix system to solve Ax=b

%Get dimensions of A
N = length(flux_data);

%Create vector of ones
col = ones(N,1);
%Create matrix A: an NxN matrix with ones along the relevant diagonals
A = spdiags([col col col col col],[-2,-1,0,1,2],N,N);

%Create vector b: give flux data a new variable name so it corresponds to
%the variable name from the matrix system in the title of this section
b = flux;

%% Solve matrix system to obtain actual beam flux values
%Use the backslash method to solve this system of equations
x = A\b;

%Create a more useful variable name
real_flux_data = x;
%% Fit Gaussian to the real flux data to obtain the Gaussian parameters
[sigma_real, mu_real] = gaussfit(beam_pos, real_flux_data);

%% Calculate the Gaussian function (making sure we scale the Gaussian correctly)

%Gaussian curve
flux_real = 1/(sqrt(2*pi)* sigma_real ) * exp( - (beam_pos-mu_real).^2 / (2*sigma_real^2));
%Compute scale factor
scale_fac = max(real_flux_data)/max(flux_real);
%Scale the gaussian
flux_real = scale_fac*flux_real;

% %Fit Gaussian using the least squares fit
% flux_real = amplitude_real * exp( - (beam_pos-mu_real).^2 / (2*sigma_real^2));

%% Calculate the FWHM for the Gaussians
FWHM_real = 2*sigma_real*sqrt(2*log(2));
FWHM = 2*sigma*sqrt(2*log(2));

%% Get correct x length for gon data
beam_pos_gon = linspace(beam_pos(1),beam_pos(end),length(gon_data(:,3)));

%% Plot all of the data in the Figures

figure('name','Flux Plots')
plot(beam_pos,flux_data,'o',beam_pos,real_flux_data,'x',beam_pos,flux,'r',beam_pos,flux_real,'k',beam_pos_gon,gon_data(:,3),'^')
title('Amp reading against aperature Position','FontSize',14)
ylabel('ipin reading (microamps)','FontSize',14)
xlabel('Aperature position (millimetres)','FontSize',14)
h = legend('Original data','Deconvoluted data','gaussian fit (Original data)','Gaussian Fit (Real Data)','Wire Scan Reading');
set(h,'FontSize',12)


end

