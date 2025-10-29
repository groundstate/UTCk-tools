% https://webtai.bipm.org/api/v0.2-beta/index.html
% curl -k --url "https://webtai.bipm.org/api/v0.2-beta/get-data.html?scale=utcr&lab=AUS&outfile=txt" > stuff


% New UTCr data every 7 days with at least 3 day latenccy
avgint = 7;
utcrlatency = 0;
utcrint = 7;
mjdStart = 58300;
nWeeks = 20;
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

sc0 = steerclock(clk,avgint,utcrint,0,mjdStart,G1,G2,nWeeks);
sc1 = steerclock(clk,avgint,utcrint,1,mjdStart,G1,G2,nWeeks);
sc2 = steerclock(clk,avgint,utcrint,2,mjdStart,G1,G2,nWeeks);
sc3 = steerclock(clk,avgint,utcrint,3,mjdStart,G1,G2,nWeeks);
sc4 = steerclock(clk,avgint,utcrint,4,mjdStart,G1,G2,nWeeks);
sc5 = steerclock(clk,avgint,utcrint,5,mjdStart,G1,G2,nWeeks);
sc6 = steerclock(clk,avgint,utcrint,6,mjdStart,G1,G2,nWeeks);
sc7 = steerclock(clk,avgint,utcrint,7,mjdStart,G1,G2,nWeeks);

figure(1);
plot(clk(:,1),clk(:,2),'.');
title('Effect of latency (constant offset)');
xlabel('MJD');
ylabel('CLK (ns)')
hold on;
plot(sc0(:,1),sc0(:,2),'.-');
plot(sc1(:,1),sc1(:,2),'.-');
plot(sc2(:,1),sc2(:,2),'.-');
plot(sc3(:,1),sc3(:,2),'.-');
plot(sc4(:,1),sc4(:,2),'.-');
plot(sc5(:,1),sc5(:,2),'.-');
plot(sc6(:,1),sc6(:,2),'.-');
plot(sc7(:,1),sc7(:,2),'.-');
hold off;
xlim([mjdStart mjdStart+(nWeeks-3)*7]);
ylim([-50 50]);


% Test data set 2
% Constant offset + noise

clk(:,2) = clk(:,2)  +(-5 + 10*rand(nWeeks*7,1));
sc0 = steerclock(clk,avgint,utcrint,0,mjdStart,G1,G2,nWeeks);
sc1 = steerclock(clk,avgint,utcrint,1,mjdStart,G1,G2,nWeeks);
sc2 = steerclock(clk,avgint,utcrint,2,mjdStart,G1,G2,nWeeks);
sc3 = steerclock(clk,avgint,utcrint,3,mjdStart,G1,G2,nWeeks);
sc4 = steerclock(clk,avgint,utcrint,4,mjdStart,G1,G2,nWeeks);
sc5 = steerclock(clk,avgint,utcrint,5,mjdStart,G1,G2,nWeeks);
sc6 = steerclock(clk,avgint,utcrint,6,mjdStart,G1,G2,nWeeks);
sc7 = steerclock(clk,avgint,utcrint,7,mjdStart,G1,G2,nWeeks);
legend('unsteered','0','1','2','3','4','5','6','7');

figure(2);
plot(clk(:,1),clk(:,2),'.');
title('Effect of latency (constant offset + noise)');
xlabel('MJD');
ylabel('CLK (ns)')
hold on;
plot(sc0(:,1),sc0(:,2),'.-');
plot(sc1(:,1),sc1(:,2),'.-');
plot(sc2(:,1),sc2(:,2),'.-');
plot(sc3(:,1),sc3(:,2),'.-');
%plot(sc4(:,1),sc4(:,2),'.-');
%plot(sc5(:,1),sc5(:,2),'.-');
%plot(sc6(:,1),sc6(:,2),'.-');
plot(sc7(:,1),sc7(:,2),'.-');
hold off;
xlim([mjdStart mjdStart+(nWeeks-3)*7]);
ylim([-50 50]);
legend('unsteered', '0','1','2','3','7');

function sc = steerclock(freeclk,avgint,utcrint,utcrlatency,mjdStart,G1,G2,nWeeks)
    sc = freeclk;
   
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
