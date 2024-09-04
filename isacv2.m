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


as = 0.9;
aw = 0.1;

h = (randn(2,2,nmc) + 1i * randn(2,2,nmc)) / sqrt(2);

w = zeros(2,2,nmc);

dr_u1 = zeros(1,nmc);
dr_u2 = zeros(1,nmc);

dr1_sorted = zeros(length(snrdbs),nmc);
dr2_sorted = zeros(length(snrdbs),nmc);

for j = 1:nmc

    w(:,1,j) = h(:,1,j)' / norm(h(:,1,j));
    w(:,2,j) = h(:,2,j)' / norm(h(:,2,j));

end

for i = 1:length(plw)

    h(:,1,:) = sqrt(pls(i)) * h(:,1,:);
    h(:,2,:) = sqrt(plw(i)) * h(:,2,:);

    for j = 1:nmc

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

dof1 = -pi:pi/8:pi;
dof2 = pi:-pi/8:-pi;

ath1 = [ones(1,length(dof1));exp(-1i*pi.*sin(dof1))];
ath2 = [ones(1,length(dof1));exp(-1i*pi.*sin(dof2))];

%covariance matrix 

Rw = w(:,1,852)'*w(:,1,852);

pth1 = zeros(1,length(ath1));
pth2 = zeros(1,length(ath2));

for i = 1:length(ath1)  
    pth1(i) = ath1(:,i)'*Rw*ath1(:,i);
    pth2(i) = ath2(:,i)'*Rw*ath2(:,i);
end


figure(1)
plot(dr1_sorted(1,:),ycdf,'g');
hold on;
plot(dr2_sorted(1,:),ycdf,'r');xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;

plot(dr1_sorted(2,:),ycdf,'g');
hold on;
plot(dr2_sorted(2,:),ycdf,'r');xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;


plot(dr1_sorted(3,:),ycdf,'g');
hold on;
plot(dr2_sorted(3,:),ycdf,'r');xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;


plot(dr1_sorted(4,:),ycdf,'g');
hold on;
plot(dr2_sorted(4,:),ycdf,'r');xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;


plot(dr1_sorted(5,:),ycdf,'g');
hold on;
plot(dr2_sorted(5,:),ycdf,'r');xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;


plot(dr1_sorted(6,:),ycdf,'g');
hold on;
plot(dr2_sorted(6,:),ycdf,'r');xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;

plot(dr1_sorted(7,:),ycdf,'g');
hold on;
plot(dr2_sorted(7,:),ycdf,'r');xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;

plot(dr1_sorted(8,:),ycdf,'g');
hold on;
plot(dr2_sorted(8,:),ycdf,'r');xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;

plot(dr1_sorted(9,:),ycdf,'g');
hold on;
plot(dr2_sorted(9,:),ycdf,'r');xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;

plot(dr1_sorted(10,:),ycdf,'g');
hold on;
plot(dr2_sorted(10,:),ycdf,'r');xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;
hold on;

plot(dr1_sorted(11,:),ycdf,'g');
hold on;
plot(dr2_sorted(11,:),ycdf,'r');xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;

figure(2)
plot(mean(dr1_sorted),ycdf,'g');
hold on;
plot(mean(dr2_sorted),ycdf,'r');xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;

