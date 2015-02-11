%% Function to calculate the rmsd of calculated values from a gaussain with experimental values using given parameters
function val_to_min = func_to_min2(noise,flux_2d)
%% State aperture diameter in microns
aperture_diameter =10;

%% Get correct kernel
kernel = get_kernel(aperture_diameter);

%% Deconvolute Signal
flux_2d_calc = deconvwnr(flux_2d,kernel,noise);

%% Find 2-norm of the difference of the matrices
val_to_min = norm(flux_2d - flux_2d_calc);

end