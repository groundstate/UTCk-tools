% https://webtai.bipm.org/api/v0.2-beta/index.html
% curl -k --url "https://webtai.bipm.org/api/v0.2-beta/get-data.html?scale=utcr&lab=AUS&outfile=txt" > stuff


% New UTCr data every 7 days
avgint = 7;
utcrlatency = 3;
utcrint = 7;
mjdStart = 58288;

clk = load('utcraus.txt');
% Start at 58288, near offset of zero so we're not trying to steer too hard
clk = clk(1240:end,:);
nWeeks=217;
G1=0.6; % gain for clock drift correction
G2=0.75;% gain for clock offset correction (via steer)

[sc,fs] = steerclock(clk,avgint,utcrint,utcrlatency,mjdStart,G1,G2,nWeeks);
sc = sc(1:end-(utcrlatency+1),:);

figure(1);
plot(clk(:,1),clk(:,2),'.');
title('UTCr');
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
[adevsc, ~, ~, tausc] = adev(sc(:,2)*1.0E-9/86400.0,rate,tau,overlap,phase,gaps);

figure(2);
loglog(taufree,adevfree,'-');
hold on;
loglog(tausc,adevsc,'-');
legend('free run','steered');
xlabel('tau (days)');
ylabel('frac. ADEV');
hold off;

[tdevfree, ~, ~, ttaufree ] = tdev(clk(:,2),rate,tau,phase);
[tdevsc, ~, ~, ttausc] = tdev(sc(:,2),rate,tau,phase);

figure(3);
loglog(ttaufree,tdevfree,'+-');
hold on;
loglog(ttausc,tdevsc,'+-');
legend('free run','steered');
xlabel('tau (days)');
ylabel('TDEV (ns)');
hold off;

% Now sweep G1 and G2 to find an optimal value
G1step = 0.01;
G1start = 0.1;
G1stop  = 1.0;
nG1steps = round((G1stop - G1start)/G1step);

G2step = 0.01;
G2start = 0.1;
G2stop  = 1.0;
nG2steps = round((G2stop - G2start)/G2step);

avg1g2 = zeros(nG1steps+1,nG2steps+1);
sdg1g2 = zeros(nG1steps+1,nG2steps+1);
deltag1g2 = zeros(nG1steps+1,nG2steps+1);
g1x = G1start + (0:nG1steps)*G1step;
g2y = G2start + (0:nG2steps)*G2step;
for iG1=0:nG1steps
    for iG2=0:nG2steps
        G1 = G1start + iG1*G1step;
        G2 = G2start + iG2*G2step;
        [sc,~] = steerclock(clk,avgint,utcrint,utcrlatency,mjdStart,G1,G2,nWeeks);
        sc = sc(1:end-(utcrlatency+1),:);
        sd = std(sc(:,2));
        av = mean(sc(:,2));
        delta = max(sc(:,2)) - min(sc(:,2));
        %fprintf('%g %g %g %g %g\n',G1,G2,av,sd,delta);
        avg1g2(iG1+1,iG2+1) = av;
        sdg1g2(iG1+1,iG2+1) = sd;
        deltag1g2(iG1+1,iG2+1) = delta;
    end
end

xt = [0 10 20 30 40 50 60 70 80 90];
xtl = {'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'};

figure(4);
contourf(deltag1g2,10);
colorbar;
xlabel('G2 ');
ylabel('G1 ');
title('Offset: max - min (ns)');
xticks(xt);
xticklabels(xtl);
yticks(xt);
yticklabels(xtl);


figure(5);
contourf(sdg1g2);
colorbar;
xlabel('G2 ');
ylabel('G1 ');
title('Offset: std dev(ns)');
xticks(xt);
xticklabels(xtl);
yticks(xt);
yticklabels(xtl);

figure(6);
contourf(avg1g2);
colorbar;
xlabel('G2 ');
ylabel('G1 ');
title('Offset: average (ns)');
xticks(xt);
xticklabels(xtl);
yticks(xt);
yticklabels(xtl);

function [sc,fs] = steerclock(freeclk,avgint,utcrint,utcrlatency,mjdStart,G1,G2,nWeeks)
    sc = freeclk;
    fs = zeros(nWeeks-2,1);
    for n=0:(nWeeks-3)
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
        fs(n+1) = fslew;
        %fprintf('%g %g %g %g %g\n', ...
        %    utcr(imjd0,1), utcr(imjd1,1),fzero,fslew,utcr(imjd1,2))
         % Now apply this to the next 7 days of data AFTER the latency offset
        % This is what we'll see in UTCr
        fsteer = G1*fzero + G2*fslew;
        for i=1:7
            % New value = current offset + underlying clock drift + steer
            sc(imjd1+i+ utcrlatency,2) = sc(imjd1+utcrlatency,2) + (freeclk(imjd1+i+utcrlatency,2) - freeclk(imjd1+utcrlatency,2)) - fsteer * i;
        end
    end
end
