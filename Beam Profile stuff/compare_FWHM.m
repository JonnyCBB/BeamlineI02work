function [FWHM_mini FWHM_apw FWHM_wire gauss_reading_mini gauss_reading_ap] = compare_FWHM(miniap_data,ap_wire_data,wire_data)

%% Obtain the beam position
beam_pos = miniap_data(:,1);
beam_pos2 = ap_wire_data(:,1);
beam_pos_gon = wire_data(:,1);

%% Fit Gaussian to the data

[sigma_mini, mu_mini] = gaussfit(beam_pos, miniap_data(:,2));
[sigma_ap, mu_ap] = gaussfit(beam_pos2, ap_wire_data(:,2));
[sigma_wire, mu_wire] = gaussfit(beam_pos_gon, wire_data(:,3));

%% Calculate the Gaussian functions (making sure we scale the Gaussian correctly)

%Gaussian curve
gauss_reading_mini = 1/(sqrt(2*pi)* sigma_mini ) * exp( - (beam_pos-mu_mini).^2 / (2*sigma_mini^2));
%Compute scale factor
scale_mini = max(miniap_data(:,2))/max(gauss_reading_mini);
%Scale the gaussian
gauss_reading_mini = scale_mini*gauss_reading_mini;

%Gaussian curve
gauss_reading_ap = 1/(sqrt(2*pi)* sigma_ap ) * exp( - (beam_pos-mu_mini).^2 / (2*sigma_ap^2));
%Compute scale factor
scale_ap = max(ap_wire_data(:,2))/max(gauss_reading_ap);
%Scale the gaussian
gauss_reading_ap = scale_ap*gauss_reading_ap;

%Gaussian curve
gauss_reading_wire = 1/(sqrt(2*pi)* sigma_wire ) * exp( - (beam_pos_gon-mu_wire).^2 / (2*sigma_wire^2));
%Compute scale factor
%scale = wire_data(round(length(wire_data(:,3))/2)-7,3)/max(gauss_reading_wire);
scale_wire = max(wire_data(:,3))/max(gauss_reading_wire);
%Scale the gaussian
gauss_reading_wire = scale_wire*gauss_reading_wire;

%% Get correct x length for gon data
beam_pos_gon = linspace(beam_pos(1),beam_pos(end),length(wire_data(:,3)));

% %% Plot the graphs
% figure('name','Gaussian fits to data')
% plot(beam_pos,gauss_reading_mini,'b',beam_pos2,gauss_reading_ap,'k',beam_pos_gon,gauss_reading_wire,'r',beam_pos,miniap_data(:,2),'o',beam_pos2,ap_wire_data(:,2),'x',beam_pos_gon,wire_data,'^')
% legend('miniap','ap/wire scan','wire scan','miniap data','ap/wire scan data','wire scan data')

%% Plot Graphs separately
figure('name','miniap data')
plot(beam_pos, miniap_data(:,2), 'o', beam_pos, gauss_reading_mini, 'b')

figure('name','ap/wire data')
plot(beam_pos2, ap_wire_data(:,2), 'x', beam_pos2, gauss_reading_ap , 'r')

figure('name','wire data')
plot(beam_pos_gon, wire_data(:,3), '^', beam_pos_gon, gauss_reading_wire , 'k')

%% Plot Graphs of miniap data and aperture w/ wire data
figure('name','miniap and ap/wire data')
plot(beam_pos,gauss_reading_mini,'g',beam_pos,gauss_reading_ap,'r')
ylabel('ipin reading')
xlabel('Position')
legend('miniap reading','ap and wire reading')

%% Calculate the FWHM
FWHM_mini = 2*sigma_mini*sqrt(2*log(2));
FWHM_apw = 2*sigma_ap*sqrt(2*log(2));
FWHM_wire = 2*sigma_ap*sqrt(2*log(2));

%% Calculate ratio of maximums between miniap and ap/wire scans
max_mini = scale_mini/(sqrt(2*pi)* sigma_mini );
max_ap = scale_ap/(sqrt(2*pi)* sigma_ap );
ratio = max_ap/max_mini;
end

