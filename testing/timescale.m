% https://webtai.bipm.org/api/v0.2-beta/index.html
% curl -k --url "https://webtai.bipm.org/api/v0.2-beta/get-data.html?scale=utcr&lab=AUS&outfile=txt" > stuff





% New UTCr data every 7 days
avgint = 7;
utcrlatency = 2;
utcrint = 7;
mjdStart = 58300;
nWeeks = 6;
G1=1.0; % gain for clock drift correction
G2=1.0; % gain for clock offset correction (via steer)

% Test data set 1
% Constant frequency offset 

foffset = -2; % ns/day
toffset = 3;
clk = zeros(nWeeks*7,2); % this is the free running clock wrt eg UTC
for i = 0:(nWeeks*7-1)
    clk(i+1,1)  = i+mjdStart;
    clk(i+1,2)  = -foffset*i + toffset;
end

sc = steerclock(clk,avgint,utcrint,mjdStart,G1,G2,nWeeks);

figure(1);
plot(clk(:,1),clk(:,2),'.');
title('clk - REF');
xlabel('MJD');
ylabel('ns')
hold on;
plot(sc(:,1),sc(:,2),'+');
hold off;

% Test data set 2
% Constant offset + noise

clk(:,2) = clk(:,2)  +(-5 + 10*rand(nWeeks*7,1));
G1 = 1.0;
G2 = 1.0;
sc = steerclock(clk,avgint,utcrint,mjdStart,G1,G2,nWeeks);

figure(2);
plot(clk(:,1),clk(:,2),'.');
title('clk - REF');
xlabel('MJD');
ylabel('ns')
hold on;
plot(sc(:,1),sc(:,2),'+');
hold off;

% Test data set 3
% Constant offset + rate
foffset = -2; % ns/day
frate   = 0.2; % ns/day^2

toffset = 3;
clk = zeros(nWeeks*7,2); % this is the free running clock wrt eg UTC
for i = 0:(nWeeks*7-1)
    clk(i+1,1)  = i+mjdStart;
    clk(i+1,2)  = -foffset*i + frate*i*i + toffset;
end

sc = steerclock(clk,avgint,utcrint,mjdStart,G1,G2,nWeeks);

figure(3);
plot(clk(:,1),clk(:,2),'.');
title('clk - REF');
xlabel('MJD');
ylabel('ns')
hold on;
plot(sc(:,1),sc(:,2),'+');
hold off;

% Doesn't work so well - but that's not surprising.

% Try some real data

clk = load('utcraus.txt');
% Start at 58287, near offset of zero so we're not trying to steer too hard
clk = clk(1239:end,:);
nWeeks=210;
G1=0.8;
G2=0.5;
sc = steerclock(clk,avgint,utcrint,mjdStart,G1,G2,nWeeks);

figure(4);
plot(clk(:,1),clk(:,2),'.');
title('clk - REF');
xlabel('MJD');
ylabel('ns')
hold on;
plot(sc(:,1),sc(:,2),'+');
hold off;

overlap = 1;
phase = 1;
gaps = 0;
rate = 1.0;
tau = (1:500);

[adevfree, ~, ~, taufree ] = adev(clk(:,2)*1.0E-9/86400.0,rate,tau,overlap,phase,gaps);
[adevsc, ~, ~, tausc] = adev(sc(1:1450,2)*1.0E-9/86400.0,rate,tau,overlap,phase,gaps);

figure(5);
loglog(taufree,adevfree,'-');
hold on;
loglog(tausc,adevsc,'-');
legend('free run','steered');
xlabel('tau (days)');
ylabel('frac. ADEV');
hold off;

[tdevfree, ~, ~, ttaufree ] = tdev(clk(:,2),rate,tau,phase);
[tdevsc, ~, ~, ttausc] = tdev(sc(1:1450,2),rate,tau,phase);

figure(6);
loglog(ttaufree,tdevfree,'+-');
hold on;
loglog(ttausc,tdevsc,'+-');
legend('free run','steered');
xlabel('tau (days)');
ylabel('TDEV (ns)');
hold off;

function sc = steerclock(freeclk,avgint,utcrint,mjdStart,G1,G2,nWeeks)
    sc = freeclk;
   
    for n=0:(nWeeks-2)
        mjd0 = mjdStart + n * utcrint;
        imjd0 = mjd0 - mjdStart + 1;
        mjd1  = mjd0 + avgint - 1; 
        imjd1 = mjd1 - mjdStart + 1;
        % Fit to the free running clock data   
        freeclkfit =  freeclk(imjd0:imjd1,:);
        [a,~,~,~] = linfit(freeclkfit(:,1),freeclkfit(:,2),0);
        % slope a(2) = dphi/dt so ffe is -a(2) and required offset to zero it is +a(2)
        fzero = a(2); 
        % calculate the slew required to bring the current offset back 
        % to zero over utcrint days

        fslew = sc(imjd1,2)/utcrint;
        
        %fprintf('%g %g %g %g %g\n', ...
        %    utcr(imjd0,1), utcr(imjd1,1),fzero,fslew,utcr(imjd1,2))
        % Now apply this to the next 7 days of data
        % This is what we'll see in UTCr
        fsteer = G1*fzero + G2*fslew;
        for i=1:7
            % New value = current offset + underlying clock drift + steer
            sc(imjd1+i,2) = sc(imjd1,2) + (freeclk(imjd1+i,2) - freeclk(imjd1,2)) - fsteer * i;
        end
    end
end
