%% Script to analyse beam data 

%% Import data from the selected files 
importfile('13763.dat');
miniap_x_data_fine = data(:,1:2);

importfile('13764.dat');
miniap_x_data_coarse = data(:,1:2);

importfile('13765.dat');
miniap_y_data_fine = data(:,1:2);

importfile('13766.dat');
miniap_y_data_coarse = data(:,1:2);

importfile('13771.dat');
gon_x_data_ap_fine = data(:,1:2);

importfile('13772.dat');
gon_x_data_ap_coarse = data(:,1:2);

importfile('13773.dat');
gon_y_data_ap_fine = data(:,1:2);

importfile('13774.dat');
gon_y_data_ap_coarse = data(:,1:2);

importfile('13781.dat');
gon_x_data_fine = data(:,1:2);

importfile('13782.dat');
gon_x_data_coarse = data(:,1:2);

importfile('13777.dat');
gon_y_data_fine = data(:,1:2);

importfile('13778.dat');
gon_y_data_coarse = data(:,1:2);

% %% Plot the x Data
% 
% figure('name','miniap x data')
% plot(miniap_x_data_fine(:,1),miniap_x_data_fine(:,2),'o',miniap_x_data_coarse(:,1),miniap_x_data_coarse(:,2),'x')
% 
% figure('name','wire scan x data')
% plot(gon_x_data_fine(:,1),gon_x_data_fine(:,2),'<',gon_x_data_coarse(:,1),gon_x_data_coarse(:,2),'>')
% 
% figure('name','aperture with wire x data')
% plot(gon_x_data_ap_fine(:,1),gon_x_data_ap_fine(:,2),'<',gon_x_data_ap_coarse(:,1),gon_x_data_ap_coarse(:,2),'>')
% 
% %% Plot the y Data
% 
% figure('name','miniap y data')
% plot(miniap_y_data_fine(:,1),miniap_y_data_fine(:,2),'o',miniap_y_data_coarse(:,1),miniap_y_data_coarse(:,2),'x')
% 
% figure('name','wire scan y data')
% plot(gon_y_data_fine(:,1),gon_y_data_fine(:,2),'<',gon_y_data_coarse(:,1),gon_y_data_coarse(:,2),'>')
% 
% figure('name','aperture with wire y data')
% plot(gon_y_data_ap_fine(:,1),gon_y_data_ap_fine(:,2),'<',gon_y_data_ap_coarse(:,1),gon_y_data_ap_coarse(:,2),'>')

%% Convert Wire reading into corresponding aperture signal
scan_gon_x_data_fine = diff(gon_x_data_fine(:,2));
scan_gon_x_data_fine(end+1) = scan_gon_x_data_fine(end);
scan_gon_x_data_fine = -(scan_gon_x_data_fine);
gon_x_data_fine = [gon_x_data_fine, scan_gon_x_data_fine];

scan_gon_x_data_coarse = diff(gon_x_data_coarse(:,2));
scan_gon_x_data_coarse(end+1) = scan_gon_x_data_coarse(end);
scan_gon_x_data_coarse = -(scan_gon_x_data_coarse);
gon_x_data_coarse = [gon_x_data_coarse, scan_gon_x_data_coarse];

scan_gon_y_data_fine = diff(gon_y_data_fine(:,2));
scan_gon_y_data_fine(end+1) = scan_gon_y_data_fine(end);
scan_gon_y_data_fine = -(scan_gon_y_data_fine);
gon_y_data_fine = [gon_y_data_fine, scan_gon_y_data_fine];

scan_gon_y_data_coarse = diff(gon_y_data_coarse(:,2));
scan_gon_y_data_coarse(end+1) = scan_gon_y_data_coarse(end);
scan_gon_y_data_coarse = -(scan_gon_y_data_coarse);
gon_y_data_coarse = [gon_y_data_coarse, scan_gon_y_data_coarse];

%% Compare miniap data with wire/ap scans
figure('name','x data')
plot(miniap_x_data_fine(:,1),miniap_x_data_fine(:,2),'o',miniap_x_data_coarse(:,1),miniap_x_data_coarse(:,2),'x',miniap_x_data_fine(:,1),gon_x_data_ap_fine(:,2),'<',miniap_x_data_coarse(:,1),gon_x_data_ap_coarse(:,2),'>')
legend('miniap fine x','miniap coarse x','wire and ap fine x','wire and ap coarse x')

figure('name','y data')
plot(miniap_y_data_fine(:,1),miniap_y_data_fine(:,2),'o',miniap_y_data_coarse(:,1),miniap_y_data_coarse(:,2),'x',miniap_y_data_fine(:,1),gon_y_data_ap_fine(:,2),'<',miniap_y_data_coarse(:,1),gon_y_data_ap_coarse(:,2),'>')
legend('miniap fine y','miniap coarse y','wire and ap fine y','wire and ap coarse y')

%% Plot wire scans
figure('name','wire scans x data')
plot(gon_x_data_fine(:,1),gon_x_data_fine(:,3),'o')

figure('name','wire scans y data')
plot(gon_y_data_fine(:,1),gon_y_data_fine(:,3),'o')

%% Compare all data and calculate FWHM

[FWHM_mini_x FWHM_apw_x FWHM_wire_x gauss_mini_x gauss_ap_x] = compare_FWHM(miniap_x_data_fine,gon_x_data_ap_fine,gon_x_data_fine);

[FWHM_mini_y FWHM_apw_y FWHM_wire_y gauss_mini_y gauss_ap_y] = compare_FWHM(miniap_y_data_fine,gon_y_data_ap_fine,gon_y_data_fine);

%% Calculate the Deconvoluted Signal - x
[deconv_reading_x beam_pos_x deconvFWHM_x] = beam_FWHM_linalg('13771.dat',gon_x_data_fine);

cum_flux = cumsum(deconv_reading_x);
cum_flux = flipud(cum_flux);
cum_meas_flux = cumsum(gon_x_data_ap_fine(:,2));
cum_meas_flux = flipud(cum_meas_flux);
cum_meas_flux2 = cumsum(miniap_x_data_fine(:,2));
cum_meas_flux2 = flipud(cum_meas_flux2);

x1 = linspace(0,100,length(cum_flux));
x2 = linspace(0,100,length(gon_x_data_fine));
x3 = linspace(0,100,length(cum_meas_flux));
x4 = linspace(0,100,length(cum_meas_flux2));
figure('name','Cumulative Plots x')
plot(x1,cum_flux,'r',x2,gon_x_data_fine(:,2),'o',x3,cum_meas_flux,'x',x4,cum_meas_flux2,'v')
legend('deconv reading','gon reading','ap and wire reading','aperture reading')

%% Calculate the Deconvoluted Signal - y
[deconv_reading_y beam_pos_y deconvFWHM_y] = beam_FWHM_linalg('13773.dat',gon_y_data_fine);

cum_flux = cumsum(deconv_reading_y);
cum_flux = flipud(cum_flux);
cum_meas_flux = cumsum(gon_y_data_ap_fine(:,2));
cum_meas_flux = flipud(cum_meas_flux);
cum_meas_flux2 = cumsum(miniap_y_data_fine(:,2));
cum_meas_flux2 = flipud(cum_meas_flux2);

x1 = linspace(0,100,length(cum_flux));
x2 = linspace(0,100,length(gon_y_data_fine));
x3 = linspace(0,100,length(cum_meas_flux));
x4 = linspace(0,100,length(cum_meas_flux2));
figure('name','Cumulative Plots y')
plot(x1,cum_flux,'r',x2,gon_y_data_fine(:,2),'o',x3,cum_meas_flux,'x',x4,cum_meas_flux2,'v')
legend('deconv reading','gon reading','ap and wire reading','aperture reading')