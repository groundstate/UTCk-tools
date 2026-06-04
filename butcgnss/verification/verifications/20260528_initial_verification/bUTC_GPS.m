%
% Note that the MATLAB library CGGTTS is available from
% 

%% Setup
BDS=1;
GAL= 2;
GLO=3;
GPS=4;
GNSS_NAMES = ["BDS","GAL","GLO","GPS"];
GNSS_CODES =  ["C","E","R","G"];
U_CAL_GNSS = [2.4,2.4,3.8,2.7]; % from Defraigne, one sigma
U_NAVMSG  = [0.2,0.1,0.1,1.3];  % from Defraigne, one sigma

startMJD = 61041;
stopMJD  = 61071; 
iGNSS = GPS;
GNSSCode = GNSS_CODES(iGNSS);
GNSSName = GNSS_NAMES(iGNSS);

rootDir = "~/src/UTCk-tools/butcgnss/verification/data";
leapSecs = 18;
GPSEpoch = 44244;
halfWin = 12.0;

%% Read CGGTTS data, filter and eyeball it
cg = CGGTTS(startMJD-1,stopMJD+1,rootDir+"/cggtts/",GNSSCode + "ZPT13",'NamingConvention','BIPM','RemoveBadTracks','yes');
totalTracks = length(cg.Tracks);
fprintf('Total tracks = %d\n',totalTracks);

% Check how many bad tracks were removed
fprintf('Bad tracks: removed %d\n',cg.BadTracks);

% Now filter
cg.Filter(cg.DSG,0,20*10); 
highDSG = totalTracks - length(cg.Tracks);
fprintf('Filter DSG: removed  %d\n',highDSG);

cg.Filter(cg.TRKL,750,780); 
shortTracks = totalTracks - highDSG - length(cg.Tracks);
fprintf('Filter TRKL: removed  %d\n',shortTracks);

cg.Filter(cg.SRSYS,-9999,9999); 
highSRSYS = totalTracks - highDSG - shortTracks - length(cg.Tracks);
fprintf('Filter SRSYS: removed  %d\n',highSRSYS);

% and plot REFSYS
figure(1);
plot(cg.Tracks(:,cg.MJD)+cg.Tracks(:,cg.STTIME)/86400 -startMJD,cg.Tracks(:,cg.REFSYS)*0.1,'.');
xlabel("MJD - " + num2str(startMJD));
ylabel('REFSYS (ns)');
title('Filtered REFSYS');

% Comment: looks nominal

%% Calculate mean REFSYS values in a window centred on UTC 0
refsys = [cg.Tracks(:,cg.MJD)+cg.Tracks(:,cg.STTIME)/86400,cg.Tracks(:,cg.REFSYS)*0.1];
% Not very efficient ...
i = 1;
avrefsys=zeros(stopMJD - startMJD + 1,2);
stdrefsys=zeros(stopMJD - startMJD + 1,2);
for mjd=startMJD:stopMJD  
    win = (refsys(:,1) >= mjd - halfWin/24.0) & (refsys(:,1) <= mjd + halfWin/24.0);
    xrefsys = refsys(win,:,:);
    avrefsys(i,:) = [mjd,mean(xrefsys(:,2))];
    stdrefsys(i,:) = [mjd,std(xrefsys(:,2))];
    i = i + 1;
end

% and plot it
figure(2);
plot(cg.Tracks(:,cg.MJD)+cg.Tracks(:,cg.STTIME)/86400 -startMJD,cg.Tracks(:,cg.REFSYS)*0.1,'.');
hold on;
plot(avrefsys(:,1)-startMJD,avrefsys(:,2),'yo-')
hold off;
xlabel(['MJD - ',num2str(startMJD)]);
ylabel('REFSYS (ns)');
title('Check on averaged REFSYS');
legend('raw REFSYS','averaged REFSYS');

% Comment: looks nominal. 

%% Load Circular T data
% Columns are
% MJD  UTC-UTC(PTB) uA uB
% Note that uB changes during the time of interest
% and I have taken the value from MJD 619039 from Circular T 457 rather
% than Circular T 456

utcutck = load('utc_ptb.txt');

% Resample this, linearly interpolating between Circular T reporting days
mjd = startMJD:stopMJD;
mjd = mjd';
interp_utck = interp1(utcutck(:,1),utcutck(:,2),mjd,'linear');
interp_utck = interp_utck';

figure(3);
plot(cg.Tracks(:,cg.MJD)+cg.Tracks(:,cg.STTIME)/86400 -startMJD,cg.Tracks(:,cg.REFSYS)*0.1,'.');
hold on;
plot(avrefsys(:,1)-startMJD,avrefsys(:,2),'yo-');
plot(utcutck(:,1)-startMJD,utcutck(:,2),'+-','LineWidth',2);
plot(mjd-startMJD,interp_utck,'o');
hold off;
xlabel("MJD - " + num2str(startMJD));
ylabel('ns');
xlim([0 35]);
title('Check on UTC-UTC(PTB)');
legend('raw REFSYS','averaged REFSYS','UTC-UTC(PTB)','Interpolated UTC-UTC(PTB)');

%% Load Time System Corrections
tsc = load(GNSSName + ".TimeSysCorr.txt");
% Compute 
i=1; % won't be too fancy - we know how to index into tsc to get data for a given MJD
dUTCGPS=zeros(stopMJD - startMJD + 1,2);
for mjd=startMJD:stopMJD
    % BDS
    % GAL
    % GPS
    % Don't have to worry about leap second logic any more ...
    N = mjd - GPSEpoch;
    Wn = floor(N/7);
    Dn = mod(N,7);
    % fprintf('%g %g %g\n',mjd,Wn,Dn);
    iTSC = i;
    deltaUTC1 = tsc(iTSC,2) + tsc(iTSC,3)*(Dn*86400 + leapSecs - tsc(iTSC,4) + 604800*(Wn - tsc(iTSC,5))); % per the ICD
    iTSC = i+1;
    deltaUTC2 = tsc(iTSC,2) + tsc(iTSC,3)*(Dn*86400 + leapSecs - tsc(iTSC,4) + 604800*(Wn - tsc(iTSC,5)));
    deltaUTC = (deltaUTC1+deltaUTC2)/2.0;
    % fprintf('%g %g %g\n',mjd,deltaUTC,deltaUTC2);
    dUTCGPS(i,:) = [mjd,-deltaUTC/1.0E-9]; % since tUTC = tGPS - dUTC
    i=i+1;
end

figure(4);
plot(cg.Tracks(:,cg.MJD)+cg.Tracks(:,cg.STTIME)/86400 -startMJD,cg.Tracks(:,cg.REFSYS)*0.1,'.');
hold on;
plot(avrefsys(:,1)-startMJD,avrefsys(:,2),'yo-');
plot(utcutck(:,1)-startMJD,utcutck(:,2),'+-','LineWidth',2);
plot(dUTCGPS(:,1)-startMJD,dUTCGPS(:,2),'*-','LineWidth',2);
hold off;
xlabel("MJD - " + num2str(startMJD));
ylabel('ns');
xlim([0 35]);
title('Check on bUTC');
legend('raw REFSYS','averaged REFSYS','UTC-UTC(PTB)','bUTC\_GNSS');

%% Compute UTC - bUTC_GNSS and UTC(PTB) - bUTC_GNSS
% UTCk - bUTC = (UTCk - GNSS) - (UTC - GNSS);
UTCk_bUTC = avrefsys(:,2) - dUTCGPS(:,2);
UTC_bUTC  = interp_utck' + UTCk_bUTC; % (UTC - bUTC) = UTC- UTCk + (UTCk - bUTC)

% Load Circular T estimates
UTC_bUTC_CirT = load('UTC_bUTC_GPS_cirt.txt');

figure(5);
plot(UTCk_bUTC,'.-');
hold on;
plot(UTC_bUTC_CirT(:,2),'.-');
plot(UTC_bUTC,'.-'); 
plot(interp_utck);
hold off;
legend('UTCk\_bUTC','UTC-bUTC (UTCk)','UTC-bUTC (CircularT)','UTC-UTCk (interpolated)');

% Uncertainties
i=1;

u_UTCk_bUTC = zeros(stopMJD - startMJD + 1,1);
u_UTC_bUTC  = zeros(stopMJD - startMJD + 1,1);

for mjd=startMJD:stopMJD
  [uACirT,uBCirT,tau] = FindNearestCirtU(utcutck,mjd);
  ttNoise = stdrefsys(i,2);
  uCirTExtrapolation = ClockExtrapolationError(uACirT,tau);
  uClockInstability = ClockInstability(tau);
  uCalGNSS = U_CAL_GNSS(iGNSS);
  uNavMsg = U_NAVMSG(iGNSS); 
  u_UTCk_bUTC(i) = sqrt(uBCirT^2 + ttNoise^2 + uCalGNSS^2 + uNavMsg^2);
  u_UTC_bUTC(i)  = sqrt(uCirTExtrapolation^2 + uClockInstability^2 + uACirT^2 + uBCirT^2 + ttNoise^2 + uCalGNSS^2 + uNavMsg^2);
  fprintf('%d %g %g [%g %g %g %g]\n',mjd,u_UTCk_bUTC(i),u_UTC_bUTC(i),ttNoise,uCirTExtrapolation,uClockInstability,uCirTExtrapolation);
  i = i + 1;
end

mjd = startMJD:stopMJD;
for i=1:length(UTCk_bUTC)
    fprintf('%5d % 10.5g +/- % 7.5g  % 10.5g +/- % 7.5g\n',mjd(i),UTCk_bUTC(i),u_UTCk_bUTC(i),UTC_bUTC(i),u_UTC_bUTC(i));
end

function [uA, uB, tau] = FindNearestCirtU(utcutck,mjd)
  uA = 0;
  uB = 0; 
  tau = 0;
  % Won't be fussy about checking boundaries
  % We know that there is a big enough data set
  for i=1:length(utcutck)-1
      if mjd >= utcutck(i,1) && mjd <=  utcutck(i+1,1) % Got it
        % Determine the minimum time interval
        tau = mjd  - utcutck(i,1);
        if (utcutck(i+1,1) - mjd) < tau
            tau = utcutck(i+1,1) - mjd;
        end
        uA = utcutck(i,3);
        uB = utcutck(i,4);
        % fprintf('%g [%g %g] %g %g %g\n',mjd,utcutck(i,1),utcutck(i+1,1),uA,uB,tau);
        return;
      end 
  end
end

function clkStab = ClockInstability(tau)
    clkStab = 2.0*sqrt(tau); % for Cs, in ns
end

function clkErr = ClockExtrapolationError(uA,tau)
    clkErr = sqrt(2.0)*uA*tau/5.0; % as per CCTF recommendation
end