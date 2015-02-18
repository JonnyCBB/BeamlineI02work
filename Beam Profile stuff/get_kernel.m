%% This function creates a kernel matrix that is used to convolve/deconvolve with our signal data
function kernel = get_kernel(aperture_diameter)

%% Calculate aperture radius from aperture diameter
aperture_radius = aperture_diameter/2;

%% Calculate from the radius how big "the matrix" should be and set all values to zeros
%The matrix consists of all surrounding points from the centre of the
%aperture that are within the area of the diameter - i.e. these matrix
%elements could potentially contribute to the signal reading.
if rem(aperture_radius,2) == 0
    kernel = zeros(aperture_radius+1);
else
    kernel = zeros(floor(aperture_radius));
end

%% Find the centre point of the matrix and set that equal to 1
%Get dimensions of the matrix (it should be a square matrix but we'll take
%both dimensions anyway)
dimensions_of_kernel = size(kernel);

%Find the indices of the centre of the matrix
index_of_centre_of_kernel = [(dimensions_of_kernel(1)+1)/2, (dimensions_of_kernel(2)+1)/2];

%Set the value at the centre of the matrix equal to 1
kernel(index_of_centre_of_kernel(1),index_of_centre_of_kernel(2)) = 1;

%% Find the indices of the points that lie within the area of the aperture

%Function bwdist calculates the euclidean distance of each element with the
%nearest non-zero value, in this case its the centre point
arbitrary_distance_matrix = bwdist(kernel);

%Since the distance between horizontally or vertically separated elements 
%in our matrix is 2microns I need to multiply the arbitrary_distance_matrix
%by 2 to get the real distance
distance_matrix = 2*arbitrary_distance_matrix;

[row_index, col_index] = find(distance_matrix <= aperture_radius);

%% Now we have the indices of all of the elements within the aperture area we set those values to 1

%State the number of points that contribute to the signal
no_of_contributing_points = length(row_index);
for contributing_point = 1 : no_of_contributing_points
    kernel(row_index(contributing_point),col_index(contributing_point)) = 1;
end

end

