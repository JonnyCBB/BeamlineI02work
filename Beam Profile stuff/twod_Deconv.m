%% Script for 2d deconvolution of beam data to obtain true FWHM
close all
clear all
clc

%% Section For User Input

%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA FILE INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Enter the file that contains the measured readings in the x-direction that
%you want to deconvolve
x_data_file_for_deconvolving = '14774.dat';

%Enter the file that contains the measured readings in the y-direction that
%you want to deconvolve
y_data_file_for_deconvolving = '14776.dat';

%Enter the file that contains the measured readings in the x-direction that
%You've measured with a new aperture to compare with predicted values
x_data_file_for_validation = '14760.dat';

%Enter the file that contains the measured readings in the y-direction that
%You've measured with a new aperture to compare with predicted values
y_data_file_for_validation = '14761.dat';

%Enter the file that contains the wire scan data in x-direction
x_data_wire_scan = '14762.dat';

%Enter the file that contains the wire scan data in y-direction
y_data_wire_scan = '14765.dat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%state the diameter of the aperture in microns
aperture_diameter = 10;

%state the aperture diameter in microns for making measurement predictions
aperture_diameter_predicted = 20;

%Take a guess for the best noise parameter. I would leave this at zero
noise_param_guess = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Obtain parameters for the 2d Gaussian by first getting parameters for 1d gaussian in x and then in y

%%%%%%%%%%%%%%%%%%%%%%%%% FOR x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read x data file to import the data 
importfile(x_data_file_for_deconvolving);

%Create relevant variables from the data file
%make variable for the postion of the beam
beam_pos_x = data(:,1);

%make variable for the flux data
flux_data_x = data(:,2);

%Perform Gaussian fit to data
[sigma_x, mu_x] = gaussfit(beam_pos_x, flux_data_x);

%%%%%%%%%%%%%%%%%%%%%%%%% FOR y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read x data file to import the data 
importfile(y_data_file_for_deconvolving);

%Create relevant variables from the data file
%make variable for the postion of the beam
beam_pos_y = data(:,1);

%make variable for the flux data
flux_data_y = data(:,2);

%Perform Gaussian fit to data
[sigma_y, mu_y] = gaussfit(beam_pos_y, flux_data_y);

%Now I have to obtain the "scaling" parameter which is going to be the
%maximum of the maximums of the two flux data sets
A = max([max(flux_data_x), max(flux_data_y)]);

%% Plot the Guassian Fits to show if the parameters are good or not 

%Calculate the fitted flux readings using the parameter values we found
%previously for a gaussian fit
flux_x = max(flux_data_x)*exp(-((beam_pos_x-mu_x).^2)/(2*sigma_x^2));
flux_y = max(flux_data_y)*exp(-((beam_pos_y-mu_y).^2)/(2*sigma_y^2));

%Plot the data for the x direction
figure('name','Gaussian fit to x data')
plot(beam_pos_x,flux_data_x,'o',beam_pos_x,flux_x,'r')
title('Figure showing measured data and the gaussian fit for the data in x direction','FontSize',14)
ylabel('ipin Reading (Arbitrary Units)','FontSize',14)
xlabel('Position (millimetres)','FontSize',14)
h1 = legend('miniap reading','Fitted Gaussian');
set(h1,'FontSize',12)

%Plot the data for the y direction
figure('name','Gaussian fit to y data')
plot(beam_pos_y,flux_data_y,'o',beam_pos_y,flux_y,'r')
title('Figure showing measured data and the gaussian fit for the data in y direction','FontSize',14)
ylabel('ipin reading (Arbitrary Units)','FontSize',14)
xlabel('Position (millimetres)','FontSize',14)
h2 = legend('miniap reading','Fitted Gaussian');
set(h2,'FontSize',12)

%% Create the 2d Gaussian
%Create a 2d grid for the position of the beam
[x,y] = meshgrid(beam_pos_x, beam_pos_y);

%Calculate the 2d profile of the beam using the equation for a 2d Gaussian
%with the calculated parameter values
z = A.*exp(-( ((x - mu_x).^2)./(2*sigma_x^2) + ((y - mu_y).^2)./(2*sigma_y^2)));

figure('name','2d Gaussian Beam Profile from measured miniap data')
surf(x,y,z,'LineStyle','none')
title('Beam Profile from measured data','FontSize',14)
ylabel('Position in y (millimetres)','FontSize',14)
xlabel('Position in x (millimetres)','FontSize',14)
zlabel('ipin Reading (Arbitrary Units)','FontSize',14)

%% Obtain correct 2d aperture size to form the kernel matrix for convolution
ker = get_kernel(aperture_diameter);

%% Find the best noise parameter
%Here we find the noise parameter such that when we convolve the
%deconvoluted data with the correct kernel we minimise the RMSD
[noise_param rmsd]= fminsearch(@(x) func_to_min2(x,z),noise_param_guess);

%% Deconvolute the signal with required aperture diamter
z_new = deconvwnr(z,ker,noise_param);

%% Plot deconvoluted signal

figure('name','Deconvolved 2d Gaussian Beam Profile')
surf(x,y,z_new,'LineStyle','none')
title('Deconvolved signal','FontSize',14)
ylabel('Position in y (millimetres)','FontSize',14)
xlabel('Position in x(millimetres)','FontSize',14)
zlabel('Predicted Deconvolved ipin reading (Arbitrary Units)','FontSize',14)

%% Obtain the parameters for the best fit of the 2d Gaussian by fitting the corresponding 1d gaussian paramters
%Obtain the maximum flux recorded from the gaussian beam along with
%the corresponding indices in the x and y directions
[max_y_value index_for_max_y_flux] = max(max(z_new,[],1));
[max_x_value index_for_max_x_flux] = max(max(z_new,[],2));
max_flux = max([max_y_value,max_x_value]);

%Find parameters for the best fit 1d gaussians
%x_data
[sigma_x_deconv, mu_x_deconv] = gaussfit(beam_pos_x, z_new(index_for_max_x_flux,:));
%y_data
[sigma_y_deconv, mu_y_deconv] = gaussfit(beam_pos_y, z_new(:,index_for_max_y_flux));

%% Plot the fitted gaussian for the deconvolved data

%Calculate the 2d profile of the beam using the equation for a 2d Gaussian
%with the 'deconv' parameter values
z_deconv = max_flux.*exp(-( ((x - mu_x_deconv).^2)./(2*sigma_x_deconv^2) + ((y - mu_y_deconv).^2)./(2*sigma_y_deconv^2)));

%Plot the data
figure('name','Gaussian Fitted Deconvolved 2d Gaussian Beam Profile')
surf(x,y,z_deconv,'LineStyle','none')
title('Gaussian Fitted Deconvolved signal','FontSize',14)
ylabel('Position in y (millimetres)','FontSize',14)
xlabel('Position in x (millimetres)','FontSize',14)
zlabel('predicted ipin reading (Arbitrary Units)','FontSize',14)

%% Print the results of the FWHM in the x and y directions

%Calculate the FWHM of the data
FWHM_orig_readings_x = 2*sigma_x*sqrt(2*log(2));
FWHM_orig_readings_y = 2*sigma_y*sqrt(2*log(2));

FWHM_deconv_x = 2*sigma_x_deconv*sqrt(2*log(2));
FWHM_deconv_y = 2*sigma_y_deconv*sqrt(2*log(2));

%Print the FWHM to the command window
fprintf('The FWHM of the ORIGINAL MINIAP READINGS in the x direction is %f\n',FWHM_orig_readings_x);
fprintf('The FWHM of the DECONVOLVED data in the x direction is %f\n\n',FWHM_deconv_x);

fprintf('The FWHM of the ORIGINAL MINIAP READINGS in the y direction is %f\n',FWHM_orig_readings_y);
fprintf('The FWHM of the DECONVOLVED data in the y direction is %f\n\n',FWHM_deconv_y);

%% Finished with deconvolution!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We've finished with the deconvolution part of the code. Now we move on to
% the predictive part of the script.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Predict the wire scan readings.

%THIS SECTION OF CODE MAY NOT BE REQUIRED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Pick the central values of the beam to compare with the wire scan data.
% [dim_y dim_x] = size(z_deconv);
% 
% central_index_x = round((dim_x+1)/2);
% central_index_y = round((dim_y+1)/2);
% 
% aperture_radius = aperture_diameter/2;
% 
% num_of_points = floor(aperture_radius/2);
% 
% indices_for_x = (central_index_x - num_of_points) : 1 : (central_index_x + num_of_points);
% indices_for_y = (central_index_y - num_of_points) : 1 : (central_index_y + num_of_points);
% 
% values_for_predicting_wire_scan_x = z_deconv(indices_for_x,:);
% values_for_predicting_wire_scan_y = z_deconv(:,indices_for_y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%To predict the reading from the wire scan we need to sum each column and
%sum each row. The sum function in Matlab will do this for us. However the
%result of this function is a vector of values that effectively represents
%the derivative of the wire scan reading. To obtain the actual predicted
%wire scan reading we'll need to perform a cumulative sum of this data.
predicted_wire_scan_readings_x_deriv = sum(z_deconv,1);
predicted_wire_scan_readings_y_deriv = sum(z_deconv,2);

%Perform cumulative sum to predict the wire scan reading
predicted_wire_scan_readings_x = cumsum(predicted_wire_scan_readings_x_deriv);
predicted_wire_scan_readings_y = cumsum(predicted_wire_scan_readings_y_deriv);

%To match with the actual wire scan readings we need to flip the values in
%the vector so we go in descending order of cumulative values
predicted_wire_scan_readings_x = fliplr(predicted_wire_scan_readings_x);
predicted_wire_scan_readings_y = flipud(predicted_wire_scan_readings_y);

%% Import the actual wire scan data 
%Import the wire scan data for both x and y directions
importfile(x_data_wire_scan);
wire_scan_data_x = data(:,1:2);

importfile(y_data_wire_scan);
wire_scan_data_y = data(:,1:2);

%Here we add a third column to the wire scan data that differentiates the
%wire scan data to obtain data that mimics the form of an aperture scan
%(i.e. the data has a Gaussian form)
sort_3rd_col_of_wire_scan_data_x = diff(wire_scan_data_x(:,2));
sort_3rd_col_of_wire_scan_data_x(end+1) = sort_3rd_col_of_wire_scan_data_x(end);
sort_3rd_col_of_wire_scan_data_x = -(sort_3rd_col_of_wire_scan_data_x);
wire_scan_data_x = [wire_scan_data_x, sort_3rd_col_of_wire_scan_data_x];

sort_3rd_col_of_wire_scan_data_y = diff(wire_scan_data_y(:,2));
sort_3rd_col_of_wire_scan_data_y(end+1) = sort_3rd_col_of_wire_scan_data_y(end);
sort_3rd_col_of_wire_scan_data_y = -(sort_3rd_col_of_wire_scan_data_y);
wire_scan_data_y = [wire_scan_data_y, sort_3rd_col_of_wire_scan_data_y];

%Because the position values on the wire scan data don't match with the
%aperture data we need to create an arbitrary position axis
x1 = linspace(0,100,length(predicted_wire_scan_readings_x));
x2 = linspace(0,100,length(wire_scan_data_x));
y1 = linspace(0,100,length(predicted_wire_scan_readings_y));
y2 = linspace(0,100,length(wire_scan_data_y));

%We now plot the predicted scan values on the same figure with arbitrary
%position units
figure('name','Comparison of Predicted Wire Scan Readings with Actual Readings - x data')
plot(x1,predicted_wire_scan_readings_x,'r',x2,wire_scan_data_x(:,2),'x')
title('Figure showing Comparison of Predicted measurements with Actual measurements for wire scan','FontSize',14)
ylabel('ipin reading','FontSize',14)
xlabel('Arbitrary Position Units','FontSize',14)
h1 = legend('Predicted Reading','Actual Readings');
set(h1,'FontSize',12)

figure('name','Comparison of Predicted Wire Scan Readings with Actual Readings - y data')
plot(y1,predicted_wire_scan_readings_y,'r',y2,wire_scan_data_y(:,2),'x')
title('Figure showing Comparison of Predicted measurements with Actual measurements for wire scan','FontSize',14)
ylabel('ipin reading','FontSize',14)
xlabel('Arbitrary Position Units','FontSize',14)
h1 = legend('Predicted Reading','Actual Readings');
set(h1,'FontSize',12)

%As an aside I want to find the FWHM of the differentiated wire scans to
%see if they match up
[sigma_wire_actual_x, mu_wire_actual_x] = gaussfit(wire_scan_data_x(:,1),wire_scan_data_x(:,3));
[sigma_wire_actual_y, mu_wire_actual_y] = gaussfit(wire_scan_data_y(:,1),wire_scan_data_y(:,3));

[sigma_wire_pred_x, mu_wire_pred_x] = gaussfit(beam_pos_x,predicted_wire_scan_readings_x_deriv);
[sigma_wire_pred_y, mu_wire_pred_y] = gaussfit(beam_pos_x,predicted_wire_scan_readings_y_deriv);

FWHM_wire_ac_x = 2*sigma_wire_actual_x*sqrt(2*log(2));
FWHM_wire_ac_y = 2*sigma_wire_actual_y*sqrt(2*log(2));

FWHM_wire_pred_x = 2*sigma_wire_pred_x*sqrt(2*log(2));
FWHM_wire_pred_y = 2*sigma_wire_pred_y*sqrt(2*log(2));

%% Convolve the calculated readings with a new aperture

%Get new kernel for the new aperture size
ker_predicted = get_kernel(aperture_diameter_predicted);

%Predict the measurements we would record with the new aperture
z_predicted = conv2(z_deconv,ker_predicted,'same');

%Plot the 2d Profile of the predicted readings
figure('name','Predicted Beam Readings')
surf(x,y,z_predicted)
title(strcat('Predicted profile of beam with aperture size of ',num2str(aperture_diameter_predicted),'microns'),'FontSize',14)
ylabel('Position in y','FontSize',14)
xlabel('Position in x','FontSize',14)
zlabel('Predicted Measurement','FontSize',14)

%% Predict the 1d data that we would expect to see with the aperture

%We need to import the data that we want to compare with.
%%%%%%%%%%%%%%%%%%%%%%%%% FOR x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read x data file to import the data 
importfile(x_data_file_for_validation);

%Create relevant variables from the data file
%make variable for the postion of the beam
beam_pos_x_for_comparison = data(:,1);

%make variable for the flux data
flux_data_x_for_comparison = data(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%% FOR y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read x data file to import the data 
importfile(y_data_file_for_validation);

%Create relevant variables from the data file
%make variable for the postion of the beam
beam_pos_y_for_comparison = data(:,1);

%make variable for the flux data
flux_data_y_for_comparison = data(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%numtoadd = abs(min(flux_data_y_for_comparison)/2);
%flux_data_y_for_comparison = numtoadd+flux_data_y_for_comparison;

%Obtain the correct row/column indices from matrix of predicted values to 
%compare with actual data
vector_of_maximums_in_x = max(z_predicted,[],2);
max_flux_x = max(flux_data_x_for_comparison);
[min_value_x index_for_x_comparison] = min(abs(vector_of_maximums_in_x - max_flux_x));
if (1-vector_of_maximums_in_x(index_for_x_comparison)/max(vector_of_maximums_in_x)) > 0.02
    fprintf('THE APERTURE WAS NOT TRANSLATED OVER THE BEAM CENTRE IN THE X-DIRECTION.\n')
end
%[max_value_x index_for_x_comparison] = max(max(z_predicted,[],2));

vector_of_maximums_in_y = max(z_predicted,[],1);
max_flux_y = max(flux_data_y_for_comparison);
[min_value_y index_for_y_comparison] = min(abs(vector_of_maximums_in_y - max_flux_y));
if (1-vector_of_maximums_in_y(index_for_y_comparison)/max(vector_of_maximums_in_y)) > 0.02
    fprintf('THE APERTURE WAS NOT TRANSLATED OVER THE BEAM CENTRE IN THE Y-DIRECTION.\n\n')
end
%[max_value_y index_for_y_comparison] = max(max(z_predicted,[],1));


%Now plot the data to see the comparison
figure('name','Comparison of predicted measurement values with actual measurement values - x data')
plot(beam_pos_x_for_comparison, z_predicted(index_for_x_comparison,:),'r',beam_pos_x_for_comparison,flux_data_x_for_comparison,'x');
title('Figure showing Comparison of Predicted measurements with Actual measurements','FontSize',14)
ylabel('ipin reading','FontSize',14)
xlabel('Position in x','FontSize',14)
h1 = legend('Predicted reading','Actual Readings');
set(h1,'FontSize',12)

figure('name','Comparison of predicted measurement values with actual measurement values - y data')
plot(beam_pos_y_for_comparison, z_predicted(:,index_for_y_comparison),'r',beam_pos_y_for_comparison,flux_data_y_for_comparison,'x');
title('Figure showing Comparison of Predicted measurements with Actual measurements','FontSize',14)
ylabel('ipin reading','FontSize',14)
xlabel('Position in y','FontSize',14)
h1 = legend('Predicted reading','Actual Readings');
set(h1,'FontSize',12)

%% Align Plots

x_pred = zeros(2,length(z_predicted(index_for_x_comparison,:)));
x_pred(1,:) = z_predicted(index_for_x_comparison,:);

%Find index for the max value in data
[max_x_flux index_for_max_x_flux] = max(flux_data_x_for_comparison);
[max_x_calc index_for_max_x_calc] = max(x_pred(1,:));

translation = mode(diff(beam_pos_x_for_comparison));

num = index_for_max_x_calc - 1;

pos_max_x = beam_pos_x_for_comparison(index_for_max_x_flux);

start_pos = pos_max_x - num*translation;

for i = 1:length(x_pred(1,:))
    x_pred(2,i) = start_pos + (i-1)*translation;
end

%% Align Plots

y_pred = zeros(length(z_predicted(index_for_y_comparison,:)),2);
y_pred(:,1) = z_predicted(:,index_for_y_comparison);

%Find index for the max value in data
[max_y_flux index_for_max_y_flux] = max(flux_data_y_for_comparison);
[max_y_calc index_for_max_y_calc] = max(y_pred(:,1));

translation = mode(diff(beam_pos_y_for_comparison));

num = index_for_max_y_calc - 1;

pos_max_y = beam_pos_y_for_comparison(index_for_max_y_flux);

start_pos = pos_max_y - num*translation;

for i = 1:length(y_pred(:,1))
    y_pred(i,2) = start_pos + (i-1)*translation;
end


%% Calculate FWHM
[sigma_pred_x, mu_pred_x] = gaussfit(beam_pos_x, z_predicted(index_for_x_comparison,:));
[sigma_pred_y, mu_pred_y] = gaussfit(beam_pos_y, z_predicted(:,index_for_y_comparison));

[sigma_comp_x, mu_comp_x] = gaussfit(beam_pos_x_for_comparison, flux_data_x_for_comparison);
[sigma_comp_y, mu_comp_y] = gaussfit(beam_pos_y_for_comparison, flux_data_y_for_comparison);

FWHM_pred_x = 2*sigma_pred_x*sqrt(2*log(2));
FWHM_pred_y = 2*sigma_pred_y*sqrt(2*log(2));

FWHM_comp_x = 2*sigma_comp_x*sqrt(2*log(2));
FWHM_comp_y = 2*sigma_comp_y*sqrt(2*log(2));

fprintf('The FWHM MINIAP READINGS in the x direction is %f\n',FWHM_comp_x);
fprintf('The FWHM of the PREDICTION data in the x direction is %f\n\n',FWHM_pred_x);

fprintf('The FWHM MINIAP READINGS in the y direction is %f\n',FWHM_comp_y);
fprintf('The FWHM of the PREDICTION data in the y direction is %f\n\n',FWHM_pred_y);

%%
flux_x_comp = max(flux_data_x_for_comparison)*exp(-((beam_pos_x_for_comparison-mu_comp_x).^2)/(2*sigma_comp_x^2));
flux_y_comp = max(flux_data_y_for_comparison)*exp(-((beam_pos_y_for_comparison-mu_comp_y).^2)/(2*sigma_comp_y^2));


%% Plot data
figure('name','Comparison of predicted measurement values with actual measurement values - x data')
plot(x_pred(2,:), x_pred(1,:),'r',beam_pos_x_for_comparison,flux_data_x_for_comparison,'x',beam_pos_x_for_comparison,flux_x_comp,'g');
title('Figure showing Comparison of Predicted measurements with Actual measurements','FontSize',14)
ylabel('ipin reading','FontSize',14)
xlabel('Position in x','FontSize',14)
h1 = legend('Predicted reading','Actual Readings','Fit');
set(h1,'FontSize',12)

figure('name','Comparison of predicted measurement values with actual measurement values - y data')
plot(y_pred(:,2), y_pred(:,1),'r',beam_pos_y_for_comparison,flux_data_y_for_comparison,'x',beam_pos_y_for_comparison,flux_y_comp,'g');
title('Figure showing Comparison of Predicted measurements with Actual measurements','FontSize',14)
ylabel('ipin reading','FontSize',14)
xlabel('Position in y','FontSize',14)
h1 = legend('Predicted reading','Actual Readings','Fit');
set(h1,'FontSize',12)