%%
% Note that the MATLAB library CGGTTS is available from
% 

%% Setup
startMJD = 61040;
stopMJD  = 61072;
GNSS = 'G';
GNSSName = 'GPS';
rootDir = '~/src/UTCk-tools/butcgnss/verification/data';
leapSecs = 18;
GPSEpoch = 44244;

%% Read CGGTTS data, filter and eyeball it
cg = CGGTTS(startMJD,stopMJD,[rootDir '/cggtts/'],[GNSS 'ZPT13'],'NamingConvention','BIPM','RemoveBadTracks','yes');
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
xlabel(['MJD - ',num2str(startMJD)]);
ylabel('REFSYS (ns)');
title('Filtered REFSYS');

% Comment: looks nominal

%% Calculate mean REFSYS values in a 2 hour window centred on UTC 0
refsys = [cg.Tracks(:,cg.MJD)+cg.Tracks(:,cg.STTIME)/86400,cg.Tracks(:,cg.REFSYS)*0.1];
% Not very efficient ...
i = 1;
avrefsys=zeros((stopMJD-1) - (startMJD +1)+1,2);
for mjd=startMJD+1:stopMJD-1
    win = (refsys(:,1) >= mjd - 1.0/24.0) & (refsys(:,1) <= mjd + 1.0/24.0);
    xrefsys = refsys(win,:,:);
    avrefsys(i,:) = [mjd,mean(xrefsys(:,2))];
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

utcptb = load('utc_ptb.txt');

figure(3);
plot(cg.Tracks(:,cg.MJD)+cg.Tracks(:,cg.STTIME)/86400 -startMJD,cg.Tracks(:,cg.REFSYS)*0.1,'.');
hold on;
plot(avrefsys(:,1)-startMJD,avrefsys(:,2),'yo-');
plot(utcptb(:,1)-startMJD,utcptb(:,2),'+-','LineWidth',2);
hold off;
xlabel(['MJD - ',num2str(startMJD)]);
ylabel('ns');
xlim([0 35]);
title('Check on UTC-UTC(PTB)');
legend('raw REFSYS','averaged REFSYS','UTC-UTC(PTB)');

%% Load Time System Corrections
tsc = load([ GNSSName '.TimeSysCorr.txt']);
% Compute 
for mjd=startMJD+1:stopMJD-1
    % BDS
    % GAL
    % GPS
    % Don't have to worry about leap seconds any more ...
    N = mjd - GPSEpoch;
    Wn = floor(N/7);
    Dn = mod(N,7);
    fprintf('%g %g %g\n',mjd,Wn,Dn);
end

%% Compute UTC - bUTC_GNSS and UTC(PTB) - bUTC_GNSS


