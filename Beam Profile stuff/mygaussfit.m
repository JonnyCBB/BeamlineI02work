%% My own Gaussian curve fit function Using a least squares fit
%finds the parameters a, b, c, data that takes the gaussian form
% y = a * exp(-((x-b)^2)/(2*c^2))
%here b is equivalent to the mean and c is equivalent to the standard
%deviation with a being the relevant scaling factor.
function [a,b,c] = mygaussfit(X,Y)

%% remove data where the amp reading is zero or negative

%find the indices of the y vector where the values are greater than zero
greater_than_zero_amp_reading = find(Y > 0);

%Use only the values of x and y where the reading is positive 
ydata = Y(greater_than_zero_amp_reading);
xdata = X(greater_than_zero_amp_reading);

%% Set up matrix for least squares fit and solve system of equations Ax = d

%make xdata an nx1 vector
xdata = xdata(:);

%Set up A matrix 
A = [xdata.^2,  xdata,  ones(size(xdata))]; 

%make ydata and nx1 vector and log the values to set up the right-hand side
d = log(ydata(:));                  

%Solve the system of equation to obtain the least-squares solution for x
x = A\d;

%% From the Solution we need to extract the parameters: a,b,c
c = sqrt(-1/(2*x(1)));
b = x(2)*(c^2);
a = exp(x(3)+((b^2)/(2*(c^2))));

end

