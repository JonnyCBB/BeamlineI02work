function pythag(a,b,diam)

c = sqrt((a^2)+(b^2));
if c <= diam/2
    fprintf('This point is inside the aperture\n')
    fprintf('Distance from centrepoint is %f\n',c)
else
    fprintf('NOT INSIDE APERTURE\n')
    fprintf('Distance from centrepoint is %f\n',c)
end
    
end

