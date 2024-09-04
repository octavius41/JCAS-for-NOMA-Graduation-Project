close all;

snrdbw = 20;
snrdbs = 10;
snras = 10^(-snrdbs/20);
snraw = 10^(-snrdbw/20); 
nmc = 10000;
sigma = 10^-12;

as = 0.1;
aw = 0.9;

h = (randn(2,2,nmc) + 1i * randn(2,2,nmc)) / sqrt(2);

h(:,1,:) = (1/sqrt(snras)) * h(:,1,:);
h(:,2,:) = (1/sqrt(snraw)) * h(:,2,:);

w = zeros(2,2,nmc);

dr_u1 = zeros(1,nmc);
dr_u2 = zeros(1,nmc);

for j = 1:nmc

        w(:,1,j) = h(:,1,j)' / norm(h(:,1,j));
        w(:,2,j) = w(:,1,j);

        x1h = h(:,1,j)'*w(:,1,j);
    
        x2h = h(:,2,j)'*w(:,1,j);
    
        dr_u1(j) = log2(1 + (as*(abs(mean(x1h))^2) / (sigma^2)));
        dr_u2(j) = log2(1 + (aw*(abs(mean(x2h))^2) / ((as*(abs(mean(x2h))^2)+sigma^2))));    

end

dr1_sorted = sort(dr_u1);
dr2_sorted = sort(dr_u2);

ycdf = (1:nmc) / nmc;

figure(1)
plot(dr1_sorted,ycdf,'g');
hold on;
plot(dr2_sorted,ycdf,'r');xlabel('bps/Hz');ylabel('probability');legend('user1','user2');title('CDF');
grid on;


