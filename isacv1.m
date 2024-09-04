close all;
% snrdb : signal to noise ratio in terms of decibel
% snra : linear value of snr
% nmc : number of monte carlo trials
% sigma : noise variance interpereted in terms of SNR 
% h : channel vectors produced for 2 antennas for 2 user for nmc trials
% w : precoding vectors derived using zero forcing
% dr_u1 : data rate calculated for user 1
% dr_u2 : data rate calculated for user 2
% x1h : absolute square of channel matrix 1 multiplied with precodings
% x2h : absolute square of channel matrix 2 multiplied with precodings
% interference from user  1 substitued by exploiting SIC

% from the NOMA paper I inferred that user 2 needs to exploit SIC since it
% has weaker channel strength and user 1 will decode directly by 
% multiplying channel with precodings. 
% (questionable safed for later meetings) 

snrdbw = 30:1:40;
snrdbs = 20:1:30;
snras = 10.^(snrdbs/10);
snraw = 10.^(snrdbw/10); 
nmc = 10000;
sigma = 10^-12;

dref = 1;
dst = 100:100:1100;
dwk = 500:100:1500;
alp = 5;

plw = (dwk ./ dref).^(-alp);
pls = (dst ./ dref).^(-alp);


as = 0.1;
aw = 0.9;

h = (randn(2,2,nmc) + 1i * randn(2,2,nmc)) / sqrt(2);

w = zeros(2,2,nmc);

dr_u1 = zeros(1,nmc);
dr_u2 = zeros(1,nmc);

dr1_sorted = zeros(length(snrdbs),nmc);
dr2_sorted = zeros(length(snrdbs),nmc);


for i = 1:length(plw)

    h(:,1,:) = sqrt(pls(i)) * h(:,1,:);
    h(:,2,:) = sqrt(plw(i)) * h(:,2,:);

    for j = 1:nmc

        w(:,1,j) = h(:,1,j)' / norm(h(:,1,j));
        w(:,2,j) = w(:,1,j);

        x1h = h(:,1,j)'*w(:,1,j);
    
        x2h = (h(:,2,j)'*w(:,1,j));
    
        dr_u1(j) = log2(1 + (as*(abs(x1h)^2) / ((aw*(abs(mean(x2h))^2) + sigma^2))));
        dr_u2(j) = log2(1 + (aw*(abs(x2h)^2) / (sigma^2)));    

    end

    dr2_sorted(i,:) = sort(dr_u1);
    dr1_sorted(i,:) = sort(dr_u2);

    h(:,1,:) = sqrt(pls(i)) / h(:,1,:);
    h(:,2,:) = sqrt(plw(i)) / h(:,2,:);

end

ycdf = (1:nmc) / nmc;

% DIRECT SENSING APPROACH

%steering vector

dof1 = -pi/2:pi/10:pi/2;
dof2 = pi/2:-pi/10:-pi/2;

ath1 = [ones(1,length(dof1));exp(-1i*pi.*sin(dof1))];
ath2 = [ones(1,length(dof1));exp(-1i*pi.*sin(dof2))];
ath3 = [1;exp(-1i*pi.*sin(pi/4))];

%covariance matrix 

Rw = w(:,1,852).*w(:,1,852)';

pth1 = zeros(1,length(ath1));
pth2 = zeros(1,length(ath2));

for i = 1:length(ath1)  
    pth1(i) = ath1(:,i)'*Rw*ath1(:,i);
    pth2(i) = ath2(:,i)'*Rw*ath2(:,i);
end

pth1 = real(pth1);
pth2 = real(pth2);

pth3 = zeros(1,nmc);

for j = 1:nmc
    
    Rw = w(:,1,j).*w(:,1,j)';
    pth3(j) = ath3'*Rw*ath3;

end

pth3 = sort(real(pth3));

cdfv = pth3 / max(pth3);
ths = linspace(0, 2, nmc);

prdetect = 1 - interp1(pth3, cdfv, ths);
prdetect(1) = 1;
prdetect(nmc) = 0;
ncdf = cumsum(prdetect) / sum(prdetect);


figure(1)
plot(dr1_sorted(1,:),ycdf,'b','LineWidth',1.5);
hold on;
plot(dr2_sorted(1,:),ycdf,'r','LineWidth',1.5);xlabel('bps/Hz');ylabel('probability');legend('user1 (strong)','user2 (weak)');title('CDF of communication data rate');
grid on;
hold on;

plot(dr1_sorted(2,:),ycdf,'b','LineWidth',1.5);
hold on;
plot(dr2_sorted(2,:),ycdf,'r','LineWidth',1.5);xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF of communication data rate');
grid on;
hold on;


plot(dr1_sorted(3,:),ycdf,'b','LineWidth',1.5);
hold on;
plot(dr2_sorted(3,:),ycdf,'r','LineWidth',1.5);xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF of communication data rate');
grid on;
hold on;


plot(dr1_sorted(4,:),ycdf,'b','LineWidth',1.5);
hold on;
plot(dr2_sorted(4,:),ycdf,'r','LineWidth',1.5);xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF of communication data rate');
grid on;
hold on;


plot(dr1_sorted(5,:),ycdf,'b','LineWidth',1.5);
hold on;
plot(dr2_sorted(5,:),ycdf,'r','LineWidth',1.5);xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;


plot(dr1_sorted(6,:),ycdf,'b','LineWidth',1.5);
hold on;
plot(dr2_sorted(6,:),ycdf,'r','LineWidth',1.5);xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;

plot(dr1_sorted(7,:),ycdf,'b','LineWidth',1.5);
hold on;
plot(dr2_sorted(7,:),ycdf,'r','LineWidth',1.5);xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;

plot(dr1_sorted(8,:),ycdf,'b','LineWidth',1.5);
hold on;
plot(dr2_sorted(8,:),ycdf,'r','LineWidth',1.5);xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;

plot(dr1_sorted(9,:),ycdf,'b','LineWidth',1.5);
hold on;
plot(dr2_sorted(9,:),ycdf,'r','LineWidth',1.5);xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;

plot(dr1_sorted(10,:),ycdf,'b','LineWidth',1.5);
hold on;
plot(dr2_sorted(10,:),ycdf,'r','LineWidth',1.5);xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;

plot(dr1_sorted(11,:),ycdf,'b','LineWidth',1.5);
hold on;
plot(dr2_sorted(11,:),ycdf,'r','LineWidth',1.5);xlabel('bps/Hz');ylabel('probability');legend('user1 (strong)','user2 (weak)');title('CDF of communication data rate');
grid on;

figure(2)
plot(mean(dr1_sorted),ycdf,'b','LineWidth',1.5);
hold on;
plot(mean(dr2_sorted),ycdf,'r','LineWidth',1.5);xlabel('bps/Hz');ylabel('Cumulative Probability');legend('user1 (strong)','user2 (weak)');title('CDF of communication data rate');
grid on;

figure(3)
plot(rad2deg(dof1),pth1,'b','LineWidth',1.5);title("Power function outputs");xlabel("Angle in degrees");ylabel("Power level");legend("P(theta)");
hold on;
%plot(rad2deg(dof2),pth2,'r');title("Power function outputs");xlabel("Angle in degrees");ylabel("Power level");legend("P(theta)");
%hold off;ü
grid on;

figure(4)
plot(pth3,ycdf,'LineWidth',1.5);title("CDF of possible received signal power at certain angle");ylabel("Cumulative Probability");xlabel("Power level");legend("F[P(theta)]")
grid on;
figure(5)
plot(ths,ncdf,'LineWidth',1.5);
xlabel('Threshold Level'); ylabel('Cumulative Probability of Detection');title('CDF of Probability of Detection for Given Thresholds');grid on;
grid on;
