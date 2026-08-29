%
% Note that the MATLAB library CGGTTS is available from
% https://github.com/groundstate/time-transfer-tools
% Version 0.1.0
% Note: needs symlog 

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
refRx  = "PT13"; % reference receiver
leapSecs = 18;
GPSEpoch = 44244;
halfWin = 12.0;

%% Read CGGTTS data from reference data lab, filter and eyeball it
cg = CGGTTS(startMJD-1,stopMJD+1,rootDir+"/cggtts/",GNSSCode + "Z" + refRx,'NamingConvention','BIPM','RemoveBadTracks','yes');
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
% figure(1);
% plot(cg.Tracks(:,cg.MJD)+cg.Tracks(:,cg.STTIME)/86400 -startMJD,cg.Tracks(:,cg.REFSYS)*0.1,'.');
% xlabel("MJD - " + num2str(startMJD));
% ylabel('REFSYS (ns)');
% title("Check on filtered REFSYS: " + GNSS_NAMES(iGNSS));

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
    % fprintf("%g %g %g\n",sum(win),xrefsys(1,2),xrefsys(end,2));
    nBef = length(xrefsys);
    % Remove outliers
    rs = xrefsys(:,2);
    rs = sort(rs);
    nChop = round(0.05*length(rs),"TieBreaker", "even"); % to match python
    filterMedian = median(rs(nChop+1:end-nChop));
    filterStd = std(rs(nChop+1:end-nChop));
    bad = abs(rs - filterMedian) >= 6*filterStd;
    rs(bad,:)=[];
    avrefsys(i,:)  = [mjd,mean(rs)];
    stdrefsys(i,:) = [mjd,std(rs)];
    %fprintf('%d %g %g->%g %g %g %g %g\n',mjd,sum(bad),nBef,length(rs),filterMedian, filterStd,avrefsys(i,2),stdrefsys(i,2));
    i = i +1;
end

% Note: checked some spot values of mean REFSYS against the output of
% butcupdate.py AOK

% and plot it
figure(1);
plot(cg.Tracks(:,cg.MJD)+cg.Tracks(:,cg.STTIME)/86400 -startMJD,cg.Tracks(:,cg.REFSYS)*0.1,'.');
hold on;
plot(avrefsys(:,1)-startMJD,avrefsys(:,2),'go-')
hold off;
xlabel(['MJD - ',num2str(startMJD)]);
ylabel('REFSYS (ns)');
title("Check on averaged REFSYS:" + GNSS_NAMES(iGNSS));
legend('filtered REFSYS','averaged REFSYS','location','southeast','box','off');

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

figure(2);
plot(cg.Tracks(:,cg.MJD)+cg.Tracks(:,cg.STTIME)/86400 -startMJD,cg.Tracks(:,cg.REFSYS)*0.1,'.');
hold on;
plot(avrefsys(:,1)-startMJD,avrefsys(:,2),'go-');
plot(utcutck(:,1)-startMJD,utcutck(:,2),'+-','LineWidth',2);
plot(mjd-startMJD,interp_utck,'ko');
hold off;

xlabel("MJD - " + num2str(startMJD));
ylabel('ns');
xlim([0 35]);
title('Check on UTC - UTC(PTB)');
legend('filtered REFSYS','averaged REFSYS','UTC - UTC(PTB)','Interpolated UTC - UTC(PTB)','location','southeast','box','off');

%% Load Time System Corrections
tsc = load(GNSSName + ".TimeSysCorr.txt");
% Compute 
i=1; % Won't be too fancy - we know how to index into tsc to get data for a given MJD
     % In particular, the first entry in the file is for startMJD - 1
dUTCGPS=zeros(stopMJD - startMJD + 1,2);
for mjd=startMJD:stopMJD
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
    % fprintf('%g %g (%g,%g)\n',mjd,deltaUTC,deltaUTC1,deltaUTC2);
    dUTCGPS(i,:) = [mjd,-deltaUTC/1.0E-9]; % since tUTC = tGPS - dUTC
    i=i+1;
end

figure(3);
plot(cg.Tracks(:,cg.MJD)+cg.Tracks(:,cg.STTIME)/86400 -startMJD,cg.Tracks(:,cg.REFSYS)*0.1,'.');
hold on;
plot(avrefsys(:,1)-startMJD,avrefsys(:,2),'go-');
plot(utcutck(:,1)-startMJD,utcutck(:,2),'+-','LineWidth',2);
plot(dUTCGPS(:,1)-startMJD,dUTCGPS(:,2),'*-','LineWidth',2);
hold off;
xlabel("MJD - " + num2str(startMJD));
ylabel('ns');
xlim([0 35]);
title('Check on bUTC');
legend('raw REFSYS','averaged REFSYS','UTC - UTC(PTB)','bUTC\_GNSS','location','southeast','box','off');

%% Compute UTC - bUTC_GNSS and UTC(PTB) - bUTC_GNSS
% UTCk - bUTC = (UTCk - GNSS) - (UTC - GNSS);
UTCk_bUTC = avrefsys(:,2) - dUTCGPS(:,2);
UTC_bUTC  = interp_utck' + UTCk_bUTC; % (UTC - bUTC) = UTC- UTCk + (UTCk - bUTC)

% Load Circular T estimates
UTC_bUTC_CirT = load('UTC_bUTC_' + GNSS_NAMES(iGNSS) + '_cirt.txt');

figure(4);
plot(UTCk_bUTC,'.-');
hold on;
plot(UTC_bUTC_CirT(:,2),'.-');
plot(UTC_bUTC,'.-'); 
%plot(interp_utck);
hold off;
xlabel("MJD - " + num2str(startMJD));
ylabel('ns');
title('Comparison with Circular T')
%legend('UTCk-bUTC','UTC-bUTC (UTCk)','UTC-bUTC (CircularT)','UTC-UTCk (interpolated)','location','southeast','box','off');
legend('UTCk - bUTC','UTC - bUTC (CircularT)','UTC - bUTC (UTCk)','location','southeast','box','off');

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
  % fprintf('%d %g %g [%g %g %g %g]\n',mjd,u_UTCk_bUTC(i),u_UTC_bUTC(i),ttNoise,uCirTExtrapolation,uClockInstability,uCirTExtrapolation);
  i = i + 1;
end

mjd = startMJD:stopMJD;
for i=1:length(UTCk_bUTC)
    % fprintf('%5d %10.5g   +/- % 7.5g  % 10.5g  +/- % 7.5g\n',mjd(i),UTCk_bUTC(i),u_UTCk_bUTC(i),UTC_bUTC(i),u_UTC_bUTC(i));
end

%% Compare with the output of butcupdate.py
verData = load(GNSSName + ".PTB.dat"); % Columns are MJD UTCk-buTC u UTC-bUTC u
% Compute fractional differences 
fprintf('\n\nFractional differences: butcupdate.py vs this script \n');
fprintf('        UTCk - bUTC      u         UTC-bUTC         u\n\n')
for i=1:length(verData)
    dUTCk_bUTC(i) = (verData(i,2) - UTCk_bUTC(i))/verData(i,2);
    du_UTCk_bUTC(i) = (verData(i,3) - u_UTCk_bUTC(i))/verData(i,3);
    dUTC_bUTC(i) = (verData(i,4) - UTC_bUTC(i))/verData(i,4);
    du_UTC_bUTC(i) = (verData(i,5) - u_UTC_bUTC(i))/verData(i,5);
    fprintf('%5d %12.3e %12.3e %12.3e %12.3e\n',verData(i,1),dUTCk_bUTC(i),du_UTCk_bUTC(i),dUTC_bUTC(i),du_UTC_bUTC(i));
end

figure(5);
plot(dUTCk_bUTC,'.'); 
hold on;
plot(du_UTCk_bUTC,'o');
plot(dUTC_bUTC,'+');
plot(du_UTC_bUTC,'*');
hold off;
symlog(gca,'y',-17);% pseudolog plot

xlabel("MJD - " + num2str(startMJD));
ylabel('fractional difference');
title('Comparison of this script with butcupdate.py');
legend('UTCk - bUTC','u(UTCk-bUTC)','UTC-bUTC','u(UTC-bUTC)','location','southeast','box','off');

% A bit redundant, but just for completeness, compare butcupdate.py with
% Circular T

fprintf('\n\nUTC - bUTC: difference between butcupdate.py and CirT (ns)\n\n');
for i=1:length(verData)
    dUTC_cirt(i) = (verData(i,4) - UTC_bUTC_CirT(i,2));
    fprintf('%5d %6.3f\n',verData(i,1),dUTC_cirt(i));
end
% and some stats
fprintf('\nmean difference = %g\n',mean(dUTC_cirt));
fprintf('std dev = %g\n',std(dUTC_cirt));

figure(6);
plot(dUTC_cirt,'+-');
xlabel("MJD - " + num2str(startMJD));
ylabel('ns');
title('Difference between butcupdate.py and Circular T');

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
