clear
%read observation file
fprintf('Reading Observation File...\n\n')
fid = fopen('site0900.01o');
if fid == -1
    errordlg(['The file does not exist.']);
    return;
end
%Read Header of Observation File
Header=30;
x=1;
while x<=Header
    [str]=fgetl(fid);
    switch x
        case x==11
            O_header(1).ApproxPositionX = str2num(str(2:14)); 
            O_header(1).ApproxPositionY = str2num(str(16:28));
            O_header(1).ApproxPositionZ = str2num(str(31:42));
        case x==12
            O_header(1).AntennaDeltaH = str2num(str(9:14)); 
            O_header(1).AntennaDeltaE = str2num(str(23:28));
            O_header(1).AntennaDeltaN = str2num(str(37:42));
        case x==13
            O_header(1).WavvelengthL1 = str2num(str(6:8)); 
            O_header(1).WavvelengthL2 = str2num(str(12:14));    
    end
    x=x+1;
end

%Check End of Header
end_of_header = 0;
while end_of_header == 0
    [str] = fgetl(fid);
    if strfind([str],'END OF HEADER')
        end_of_header=1;
    end
end

% Read Content of Observations file
j=1;
 while feof(fid) ~= 1 %Loop unitl the file end
%Skip the headers stored in the middle of observation file
 flag=0;
 while flag==0
     [str]=fgetl(fid);
     test=str2num(str(1:4));
     a=isempty(test);
     b=strfind(str,'COMMENT');
     if a | b
         flag=0;
     else
         flag=1;
     end
 end
 
 %Store the first line of epoch
 O_eph(j).year = str2num(str(2:3)); 
 O_eph(j).month = str2num(str(5:6)); 
 O_eph(j).day = str2num(str(8:9)); 
 O_eph(j).hour = str2num(str(11:12)); 
 O_eph(j).min = str2num(str(14:15)); 
 O_eph(j).sec = str2num(str(17:27)); 
 
 %Convert calender date to GPS time using fuction
 O_eph(j).GPStime= ConvertCalenderDateToGPSTime(O_eph(j).year,O_eph(j).month,O_eph(j).day,O_eph(j).hour,O_eph(j).min,O_eph(j).sec);

 %Read number of satellite in each epoch
 O_eph(j).NumberOfSatellite=str2num(str(30:32));
 NumberOfSatellite=O_eph(j).NumberOfSatellite;
 
 %Read Satellite PRN
 x=1;
 start=34;
 ending=35;
 while x<=NumberOfSatellite
     O_eph(j).sat(x)=str2num(str(start:ending));
     start=start+3;
     ending=ending+3;
     x=x+1;
 end
 
 %Read Observations
 y=0;
while y<NumberOfSatellite
    [str]=fgetl(fid);
    start2=1;
    ending2=14;
    x=1;
    while x<=5
        O_eph(j).obs(y*5+x)=str2num(str(start2:ending2));
        start2=start2+16;
        ending2=ending2+16;
        x=x+1;
    end 
y=y+1;
[str]=fgetl(fid);
end
j=j+1;
 end
 
fprintf('The Observation File Imported\n\n')

%
%%read NAV
% PRN    ....... satellite PRN          
% M0     ....... mean anomaly at reference time
% delta_n  ..... mean motion difference
% e      ....... eccentricity
% sqrt(A)  ..... where A is semimajor axis
% OMEGA  ....... LoAN at weekly epoch
% i0     ....... inclination at reference time
% omega  ....... argument of perigee
% OMEGA_dot  ... rate of right ascension 
% i_dot  ....... rate of inclination angle
% Cuc    ....... cosine term, arg. of latitude
% Cus    ....... sine term, arg. of latitude
% Crc    ....... cosine term, radius
% Crs    ....... sine term, radius
% Cic    ....... cosine term, inclination
% Cis    ....... sine term, inclination
% toe    ....... time of ephemeris
% IODE   ....... Issue of Data Ephemeris
% GPS_wk ....... GPS week

fid = fopen('site0900.01n');
if fid == -1
    errordlg(['The file does not exist.']);
    return;
end

fprintf('Reading the Navigation File...\n\n')

%Read Header Of Navigation File
Header=7;
x=1;
while x<=Header
    [str]=fgetl(fid);
    switch x
        case x==4
        N_header(1).A0 = str2num(str(4:14)); 
        N_header(1).A1 = str2num(str(16:26));
        N_header(1).A2 = str2num(str(27:38));
        N_header(1).A3 = str2num(str(39:50));
        
        case x==5 %ION BETA
        N_header(1).B0 = str2num(str(4:14)); 
        N_header(1).B1 = str2num(str(16:26));
        N_header(1).B2 = str2num(str(27:38));
        N_header(1).B3 = str2num(str(39:50));
        
        case x==6 % Almanac parameters to compute time in UTC
        N_header(1).UTC_AO = str2num(str(4:22));  % A0: terms of polynomial  
        N_header(1).UTC_AI = str2num(str(23:41)); % A1: terms of polynomial  
        N_header(1).UTC_T = str2num(str(46:50)); %T : reference time for UTC data  
        N_header(1).UTC_W = str2num(str(56:59)); % W: UTC reference week number 
    end
    x=x+1;
end

%Check End of Header
end_of_header = 0;
while end_of_header == 0
    current_line = fgetl(fid);
    if strfind(current_line,'END OF HEADER')
        end_of_header=1;
    end
end

%Read Content of Navigation File
j = 1;
while feof(fid) ~= 1
% first line of each epoch
[str]=fgetl(fid);
N_eph(j).PRN = str2num(str(1:2)); 
N_eph(j).year=str2num(str(4:5));
N_eph(j).month=str2num(str(7:8));
N_eph(j).day=str2num(str(10:11));
N_eph(j).hour=str2num(str(13:14));
N_eph(j).minute=str2num( str(16:17));
N_eph(j).second=str2num(str(18:22));
%Convert calender date to GPS time using fuction
N_eph(j).GPStime= ConvertCalenderDateToGPSTime(N_eph(j).year,N_eph(j).month,N_eph(j).day,N_eph(j).hour,N_eph(j).minute,N_eph(j).second);  

N_eph(j).a_f0 = str2num(str(23:41));
N_eph(j).a_f1 = str2num(str(42:60));
N_eph(j).a_f2 = str2num(str(61:79));  
 
% 2
[str]=fgetl(fid);
N_eph(j).IODE_sf2 = str2num(str(4:22));   
N_eph(j).C_rs = str2num(str(23:41));
N_eph(j).deltan = str2num(str(42:60));   
N_eph(j).M_0 = str2num(str(61:79));

% 3 
[str]=fgetl(fid);
N_eph(j).C_uc = str2num(str(4:22));
N_eph(j).e = str2num(str(23:41));
N_eph(j).C_us = str2num(str(42:60)); 
N_eph(j).sqrtA = str2num(str(61:79));

% 4 
[str]=fgetl(fid);
N_eph(j).toe = str2num(str(4:22));
N_eph(j).cic = str2num(str(23:41));
N_eph(j).omega_0 = str2num(str(42:60)); 
N_eph(j).cis = str2num(str(61:79));

% 5   
[str]=fgetl(fid);
N_eph(j).io = str2num(str(4:22));
N_eph(j).Crc = str2num(str(23:41));
N_eph(j).omega = str2num(str(42:60)); 
N_eph(j).omegaDot = str2num(str(61:79));

%6 
[str]=fgetl(fid);
N_eph(j).iDot = str2num(str(4:22));
N_eph(j).weekNumber = str2num(str(42:60));  

% 7       
[str]=fgetl(fid);
N_eph(j).accracy = str2num(str(4:22));
N_eph(j).health = str2num(str(23:41));
N_eph(j).T_GD = str2num(str(42:60)); 
N_eph(j).IODC = str2num(str(61:79));

 % 8
[str]=fgetl(fid);
 N_eph(j).tST = str2num(str(4:22)); %t, gps signal transmission time in unit of seconds
j = j+1;
end

fprintf('The Navgation File Imported.\n\n')

fprintf('Matching Navigation Message to Observation File...\n\n')

%Match NAV to OBS
 j=1;
while j<=length(O_eph)
checkTIME= O_eph(j).GPStime;
satnumber=O_eph(j).NumberOfSatellite;
z=1;
while z<=satnumber
    checkPRN= O_eph(j).sat(z);
    x=1;
    flag=0;
    temptime=7200;
    clear tempx;
    %Check is each navigation match the GPS time(within four hours) 
    %and PRN of each observations 
    while x<=length(N_eph)
        if checkPRN==N_eph(x).PRN
         temp=N_eph(x).GPStime;
         test=temp-checkTIME; 
         %Test=The time difference between navigation message and observation
         if test<=14400 & test<=temptime & test>=0 
         %test<=temptime help find the navigation message with smallest
         %time difference(most updated navigation message)
             temptime=test;
             tempx=x;
             flag=1;
         end
        end
         x=x+1;
    end
    %When flag=0, no navigation message is suitable to the observation.
    if flag==0
        O_eph(j).navNO(z)=0;
    end
    %When flag=1, the position of navigation message suitable will be
    %stored
    if flag==1
        starting=3;
        if z==1
            O_eph(j).navNO(z)=tempx;
        else
            a=starting+((z-1)*5);
            O_eph(j).navNO(z)=tempx;
        end
    end
    z=z+1;
end
j=j+1;    
end
fprintf('Files Matched.\n\n')

fprintf('Combining the Navigation file and the Observation file...\n\n')
%Create the structure array COM which combines the Navigation File and Observation File 
j=1;
count=0;%Count the cuurect amount of observation that have matched with navigation message
while j<=length(O_eph)
    satnumber=O_eph(j).NumberOfSatellite;
    z=1;
    while z<=satnumber
        count=count+1;
        n=O_eph(j).navNO(z);
        if ~n==0
            %Store the observation data into Com
            Com(count).PRN=O_eph(j).sat(z);
            Com(count).year=O_eph(j).year;
            Com(count).month=O_eph(j).month;
            Com(count).day=O_eph(j).day;
            Com(count).hour=O_eph(j).hour;
            Com(count).min=O_eph(j).min;
            Com(count).sec=O_eph(j).sec;
            Com(count).GPStime=O_eph(j).GPStime;
            %Store the Corresponding Navagation Data into Com
            Com(count).NAVNO=n;
            Com(count).NAVGPStime=N_eph(n).GPStime;
            Com(count).a_f0 =N_eph(n).a_f0; 
            Com(count).a_f1  =N_eph(n).a_f1 ;
            Com(count).a_f2 =N_eph(n).a_f2 ;
            Com(count).IODE_sf2 =N_eph(n).IODE_sf2;   
            Com(count).C_rs  =N_eph(n).C_rs ;
            Com(count).deltan =N_eph(n).deltan ;
            Com(count).M_0  =N_eph(n).M_0; 
            Com(count).C_uc  =N_eph(n).C_uc ;
            Com(count).e =N_eph(n).e;
            Com(count).C_us =N_eph(n).C_us;
            Com(count).sqrtA  =N_eph(n).sqrtA ;
            Com(count).toe  =N_eph(n).toe ;
            Com(count).cic =N_eph(n).cic;
            Com(count).omega_0 =N_eph(n).omega_0 ; 
            Com(count).cis  =N_eph(n).cis ;
            Com(count).io =N_eph(n).io ;
            Com(count).Crc  =N_eph(n).Crc; 
            Com(count).omega =N_eph(n).omega; 
            Com(count).omegaDot =N_eph(n).omegaDot;
            Com(count).iDot =N_eph(n).iDot;
            Com(count).weekNumber =N_eph(n).weekNumber;
            Com(count).accracy =N_eph(n).accracy;
            Com(count).health  =N_eph(n).health ;
            Com(count).T_GD =N_eph(n).T_GD ;
            Com(count).IODC =N_eph(n).IODC;
            Com(count).tST  =N_eph(n).tST ;
            %Store the observation into the Com
            starting=3;
            if z==1
                Com(count).obs=O_eph(j).obs(starting);
            else
                Com(count).obs=O_eph(j).obs(starting+(z-1)*5);
            end
        else
                    count=count-1;
        end
    z=z+1;           
    end
j=j+1;
end
fprintf('Navigation file and Observation file Combined.\n\n')

fprintf('Caculating the Satellite Position...\n\n')

%Calculate the Satellite Position
num = 1;
while num <=length(Com)
%Define the constant
u = 3.98600418e14; % Gravitational constant (GM, in %m^3/s^2)
c= 2.99792458e8;%Speed of light (m/s)    
We=7.2921151467e-5; %earth rotaton rate

%Step1. Calculate the mean motion n
%1.	Compute the semi-major axis 
a = (Com(num).sqrtA)^2; 

%Get from the Navigation File
toc=Com(num).NAVGPStime;
toe = (Com(num).toe); 
af0=Com(num).a_f0;
af1=Com(num).a_f1;
af2=Com(num).a_f2;
delta_n = Com(num).deltan;
M0 = (Com(num).M_0); 
e = (Com(num).e); 
omega = (Com(num).omega); 
C_us = Com(num).C_us; 
C_uc= Com(num).C_uc;
C_rs = Com(num).C_rs;
C_rc = Com(num).Crc;

%Get from the Observation file
ti= Com(num).GPStime;
Rpi=Com(num).obs;

% 3. Signal transmission time 
t=ti-(Rpi/c);
dt=t-toc;
saterror=af0+af1*dt+af2*(dt^2);
Com(num).saterror=saterror;
t=t-Com(num).saterror;
Com(num).t=t;

%4.compute the tk, i.e. the GPS orbit epoch t 
%with respect to the time of ephemerides (toe)
tk = t - toe;
while tk>302400
    tk=tk-604800;
end
while tk<-302400
    tk=tk+604800;
end

% mean motion
n0 = sqrt(u/(a)^3); %(in radian/second) , or sqrt(GM/(a^3))
n=n0+delta_n;

%Step2. Calculate the mean anomaly M 
M =M0+ n * tk;
% to ensure it is within the 0 degree - 360 degrees
M = rem (M+2*pi,2*pi); %in radian

%Step 3.Calculate the eccentric anomaly E
% Initial guess of eccentric anomaly
E=M;
    for ii = 1:10
     E= M + e * sin(E);
    end
    
% to ensure it is within the 0 degree - 360 degrees
E= rem(E+2*pi,2*pi);

%Step4. Calculate the true anomaly £c (in radian)
theta = atan2 ( (sin(E)*sqrt(1-e^2)), (cos(E)-e) ); 

% Step5.Calculate the orbit radius r
r=a*(1-e*cos(E));

phi = theta + omega; % Argument of latitude in radian
phi = rem (phi+2*pi, 2*pi); % To ensure it is within the 0 degree - 360 degrees

delta_phi = C_us * sin(2*phi) + C_uc * cos(2*phi);
phi = phi + delta_phi;

delta_r = C_rs*sin(2*phi)+C_rc*cos(2*phi);
r=r+delta_r ;

%Calculate the GPS satellite position (x, y, z) in the orbital plane
%Get from Navigation File
i0 = Com(num).io; % Inclination Angle at Reference Time
IDOT = Com(num).iDot; % Rate of Inclination Angle
C_is = Com(num).cis;
C_ic = Com(num).cic;

i = i0 + IDOT * tk;% GPS Orbit Plane's Inclination
delta_i = C_is * sin(2*phi) + C_ic * cos (2*phi);
i = i + delta_i;

Com(num).ECI1=[r*cos(phi);r*sin(phi);0];

%%Satellite position Transformation
rotation=[1 0 0; 0 cos(i) -sin(i);0 sin(i) cos(i)];
ECI2=rotation*Com(num).ECI1;
Com(num).ECI2=ECI2;

asc= Com(num).omega_0+(Com(num).omegaDot- We)* tk-We*toe;
rotation2=[cos(asc) -sin(asc) 0;sin(asc) cos(asc) 0;0 0 1];
Com(num).ECEF= rotation2*Com(num).ECI2;

%sagnac effect
tt=Rpi/c-Com(num).saterror;%travellingtime
R=[cos(We*tt) sin(We*tt) 0;-sin(We*tt) cos(We*tt) 0;0 0 1];
Com(num).FinalSP=R*Com(num).ECEF;

num=num+1;
end
fprintf('The Satellite Position with correction in ECEF calculated.\n\n')

fprintf('Calculating the recevier position using Least Square Adjustment...\n\n')

% %GPS Receiver position using SPP
rec=[1;1;1;1];
for iteration=1:10
X0=rec(1);
Y0=rec(2);
Z0=rec(3);
ct=rec(4);
num = 1;
B=[];
f=[];
while num <= length(Com)
    XP=Com(num).FinalSP(1);
    YP=Com(num).FinalSP(2);
    ZP=Com(num).FinalSP(3);
    Rpi=Com(num).obs;
    saterror=Com(num).saterror;
    rho=sqrt((X0-XP)^2+(Y0-YP)^2+(Z0-ZP)^2);
    B=[B;(X0-XP)/rho (Y0-YP)/rho (Z0-ZP)/rho 1];
    f=[f;Rpi-rho+saterror*c];
    num=num+1;
end
P=eye(length(Com));
dx=[inv(B'*P*B)*B'*P*f];
rec=rec+dx;
end
rec2=[];
rec2=[rec(1) rec(2) rec(3)];
lla = ecef2lla(rec2);
fprintf('The recevier position in ECEF(m) is calculated \n\n')

fprintf('X= %.7f\n',rec(1));
fprintf('Y= %.7f\n',rec(2));
fprintf('Z= %.7f\n\n',rec(3));

fprintf('The recevier position in WGS84(d,d,m) is calculated \n\n')

fprintf('X= %.7f\n',lla(1));
fprintf('Y= %.7f\n',lla(2));
fprintf('Z= %.7f\n',lla(3));

% function GPStime=ConvertCalenderDateToGPSTime(year,month,day,hour,min,sec)
% year=year;
% month=month;
% day=day;
% hour=hour;
% min=min;
% sec=sec;
% 
% % Seconds in one week
% secs_per_week = 604800;
% 
% % Converts the two digit year to a four digit year.
% % Two digit year represents a year in the range 1980-2079.
% if (year >= 80 & year <= 99)
%     year = 1900 + year;
% end
% if (year >= 0 & year <= 79)
%     year = 2000 + year;
% end
% 
% % Calculates the 'm' term used below from the given calendar month.
% if (month <= 2)
%     y = year - 1;
%     m = month + 12;
% end
% if (month > 2)
%     y = year;
%     m = month;
% end
% 
% % Computes the Julian date corresponding to the given calendar date.
% JD = floor( (365.25 * y) ) + floor( (30.6001 * (m+1)) ) + day + ( (hour + min / 60 + sec / 3600) / 24 ) + 1720981.5;
% % Computes the GPS week corresponding to the given calendar date.
% gps_week = floor( (JD - 2444244.5) / 7 );
% % Computes the GPS seconds corresponding to the given calendar date.
%  GPStime=round(((((JD-2444244.5)/7)-gps_week)*secs_per_week)/0.5)*0.5;
% end
