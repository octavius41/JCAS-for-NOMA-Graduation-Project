close all;

% Parameters
snrdbw = 30:1:40;
snrdbs = 20:1:30;
snras = 10.^(snrdbs/10);
snraw = 10.^(snrdbw/10); 
nmc = 100;
sigma = 10^-12;
dref = 1;
dst = 100:100:1100;
dwk = 500:100:1500;
alp = 5;
plw = (dwk ./ dref).^(-alp);
pls = (dst ./ dref).^(-alp);
as = 0.9;
aw = 0.1;

% Channel and beamforming initialization
h = (randn(2,2,nmc) + 1i * randn(2,2,nmc)) / sqrt(2);
w = zeros(2,2,nmc);

% Data rate initialization
dr_u1 = zeros(1,nmc);
dr_u2 = zeros(1,nmc);
dr1_sorted = zeros(length(snrdbs),nmc);
dr2_sorted = zeros(length(snrdbs),nmc);

% Optimization parameters
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

for i = 1:length(plw)

    h(:,1,:) = sqrt(pls(i)) * h(:,1,:);
    h(:,2,:) = sqrt(plw(i)) * h(:,2,:);

    for j = 1:nmc

        % Define the objective function for optimization
        objective = @(w) -objective_function(w, h(:,:,j), as, aw, sigma);
        
        % Initial beamforming vectors
        w_init = randn(4,1);
        
        % Constraints: Power normalization and other constraints
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        lb = [];
        ub = [];
        
        % Optimize the beamforming vectors
        [w_opt, ~] = fmincon(objective, w_init, A, b, Aeq, beq, lb, ub, [], options);
        
        % Reshape w_opt back to 2x2 matrix
        w(:,:,j) = reshape(w_opt, 2, 2);
        
        % Calculate data rates for the optimized beamforming vectors
        x1h = h(:,1,j)' * w(:,1,j);
        x2h = h(:,2,j)' * w(:,1,j);
        
        dr_u1(j) = log2(1 + (as*(abs(x1h)^2) / ((aw*(abs(mean(x2h))^2) + sigma^2))));
        dr_u2(j) = log2(1 + (aw*(abs(x2h)^2) / (sigma^2)));
    end

    dr2_sorted(i,:) = sort(dr_u1);
    dr1_sorted(i,:) = sort(dr_u2);

    h(:,1,:) = sqrt(pls(i)) / h(:,1,:);
    h(:,2,:) = sqrt(plw(i)) / h(:,2,:);

end

ycdf = (1:nmc) / nmc;

% Sensing Approach
dof1 = 0:pi/10:pi;
dof2 = -pi:pi/10:0;
ath1 = [ones(1,length(dof1)); exp(-1i*pi.*sin(dof1))];
ath2 = [ones(1,length(dof2)); exp(-1i*pi.*sin(dof2))];
ath3 = [1; exp(-1i*pi.*sin(pi/4))];

pth1 = zeros(1,length(ath1));
pth2 = zeros(1,length(ath2));

for i = 1:length(ath1)
    Rw = w(:,1,85) * w(:,1,85)';
    pth1(i) = ath1(:,i)' * Rw * ath1(:,i);
    pth2(i) = ath2(:,i)' * Rw * ath2(:,i);
end

pth1 = real(pth1);
pth2 = real(pth2);

pth3 = zeros(1,nmc);

for j = 1:nmc
    Rw = w(:,1,j) * w(:,1,j)';
    pth3(j) = ath3' * Rw * ath3;
end

pth3 = sort(real(pth3));
cdfv = pth3 / max(pth3);
ths = linspace(0, 2, nmc);
prdetect = 1 - interp1(pth3, cdfv, ths);
prdetect(1) = 1;
prdetect(nmc) = 0;
ncdf = cumsum(prdetect) / sum(prdetect);

% Plotting results
figure(1)
plot(dr1_sorted(1,:), ycdf, 'g');
hold on;
plot(dr2_sorted(1,:), ycdf, 'r'); xlabel('bps/Hz'); ylabel('Probability'); legend('user1', 'user2'); title('CDF');
grid on;

for i = 2:length(plw)
    plot(dr1_sorted(i,:), ycdf, 'g');
    hold on;
    plot(dr2_sorted(i,:), ycdf, 'r'); xlabel('bps/Hz'); ylabel('Cumulative Probability'); legend('user1', 'user2'); title('CDF');
    grid on;
end

figure(2)
plot(mean(dr1_sorted), ycdf, 'g');
hold on;
plot(mean(dr2_sorted), ycdf, 'r'); xlabel('bps/Hz'); ylabel('Cumulative Probability'); legend('user1', 'user2'); title('CDF');
grid on;

figure(3)
plot(rad2deg(dof1), pth1, 'g');
hold on;
plot(rad2deg(dof2), pth2, 'r'); title("Power function outputs"); xlabel("Angle in degrees"); ylabel("Power level"); legend("P(theta)");
hold off;

figure(4)
plot(pth3, ycdf); title("CDF of possible received signal power at certain angle"); ylabel("Cumulative Probability"); xlabel("Power level"); legend("F[P(theta)]");

figure(5)
plot(ths, ncdf);
xlabel('Threshold Level'); ylabel('Cumulative Probability of Detection'); title('CDF of Probability of Detection for Given Thresholds'); grid on;

% Objective function definition
function f = objective_function(w, h, as, aw, sigma)
    w = reshape(w, 2, 2); % Reshape w back to 2x2 matrix
    x1h = h(:,1)' * w(:,1);
    x2h = h(:,2)' * w(:,1);
    dr_u1 = log2(1 + (as*(abs(x1h)^2) / ((aw*(abs(mean(x2h))^2) + sigma^2))));
    dr_u2 = log2(1 + (aw*(abs(x2h)^2) / (sigma^2)));
    sensing_power = sum(abs(w(:)).^2); % Example sensing power calculation
    f = -(dr_u1 + dr_u2 + sensing_power); % Maximize the sum of data rates and sensing power
end
