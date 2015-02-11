%% Function to calculate the Root Mean Squared Deviation (RMSD)
function rmsd_val = rmsd(y1,y2)
rmsd_val = sqrt(sum((y1-y2).^2)/length(y1-y2));
end

