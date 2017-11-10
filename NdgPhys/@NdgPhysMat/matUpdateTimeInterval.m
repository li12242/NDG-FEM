%> @brief Update the time interval
%> The default update function returns the const time interval in options
function dt = matUpdateTimeInterval( obj, fphys )
dt = obj.getOption('timeInterval');
end% func