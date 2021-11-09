function GPStime=ConvertCalenderDateToGPSTime(year,month,day,hour,min,sec)
year=year;
month=month;
day=day;
hour=hour;
min=min;
sec=sec;

% One week in seconds
secs_per_week = 604800;
%The GPS year range from 1980-2079
if (year >= 80 & year <= 99)
    year = 1900 + year;
end
if (year >= 0 & year <= 79)
    year = 2000 + year;
end
% Calculates the 'm' term used below from the given calendar month.
if (month <= 2)
    y = year - 1;
    m = month + 12;
end
if (month > 2)
    y = year;
    m = month;
end
% Computes the Julian date corresponding to the given calendar date.
JD = floor( (365.25 * y) ) + floor( (30.6001 * (m+1)) ) + day + ( (hour + min / 60 + sec / 3600) / 24 ) + 1720981.5;
% Computes the GPS week corresponding to the given calendar date.
gps_week = floor( (JD - 2444244.5) / 7 );
% Computes the GPS seconds corresponding to the given calendar date.
 GPStime=round(((((JD-2444244.5)/7)-gps_week)*secs_per_week)/0.5)*0.5;
end
