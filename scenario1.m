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

snrdbw = 40:1:50;
snrdbs = 30:1:40;
snras = 10.^(snrdbs/10);
snraw = 10.^(snrdbw/10); 
nmc = 10000;
sigma = 10^-12;

dref = 1;
dst = 100:100:1100;
dwk = 500:100:1500;
alp = 4;

plw = (dwk ./ dref).^(-alp);
pls = (dst ./ dref).^(-alp);

as = 0.9;
aw = 0.1;

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
    
        dr_u1(j) = log2(1 + (aw*(abs(mean(x1h))^2) / (sigma^2)));
        dr_u2(j) = log2(1 + (as*(abs(mean(x2h))^2) / (aw*(abs(mean(x1h))^2)  + sigma^2)));    

    end

    dr1_sorted(i,:) = sort(dr_u1);
    dr2_sorted(i,:) = sort(dr_u2);

    h(:,1,:) = sqrt(pls(i)) / h(:,1,:);
    h(:,2,:) = sqrt(plw(i)) / h(:,2,:);

end



ycdf = (1:nmc) / nmc;

figure(1)
plot(dr1_sorted(1,:),ycdf,'b','LineWidth',1.5);
hold on;
plot(dr2_sorted(1,:),ycdf,'r','LineWidth',1.5);xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF of communication data rate');
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