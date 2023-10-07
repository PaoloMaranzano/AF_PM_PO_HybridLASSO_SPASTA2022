function degrees=degreesMinutesToDecimal(degrees_minutes)
% Convert latitude and longitude between
% decimal degrees, minutes, and seconds and degrees
% Convert degrees minutes to degrees decimal
% Input is an array: [degrees,minutes,seconds] 
% 
% Example: To convert -130 40'2.5" 
%   degrees_minutes = [-130,40,2.5]
%   degrees = degreesMinutesToDecimal([130,40,2])
%   degrees =
%       130.6672
%
% Matthew Perrett 2012

if degrees_minutes(1) >= 0
    degrees = degrees_minutes(1) + degrees_minutes(2)/60 + degrees_minutes(3)/3600;
else
    degrees = degrees_minutes(1) - degrees_minutes(2)/60 - degrees_minutes(3)/3600;
end
