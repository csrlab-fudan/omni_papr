clear;
antRow = 16;
antCol = 16;
paprs = db2pow(linspace(0, 6, 31));
extendFactor = 1.025; % the circle extending factor to diminish the classification error
overSamplingRate = 10;
N1 = antRow; M1= overSamplingRate*N1;
N2 = antCol; M2 = overSamplingRate*N2;

epsilon = 1e-8;
lp = 1; % lp norm
avgs = zeros(3, length(paprs));
stds = zeros(3, length(paprs));
domains = {'circle', 'square'};
%%
for i = 1:length(domains)
    domain = domains{i};
    for j = 1:length(paprs)
        papr = paprs(j);
        X = pgd_acc(N1, N2, overSamplingRate, papr, lp, epsilon, [], domain);
        [avgs(i, j), stds(i, j)] = comp_std(X);
        if j>1 && stds(i, j)>stds(i, j-1)
            X = pgd_acc(N1, N2, overSamplingRate, papr, lp, epsilon, X_old, domain);
            [avgs(i, j), stds(i, j)] = comp_std(X);
        end
        X_old = X;
    end
end

xis = linspace(0.002, 0.4, 31);
paprs3 = zeros(size(xis));
for i = 1:length(xis)
    xi = xis(i);
    [x1, ~] = gen_qiao(xi, N1);
    [x2, ~] = gen_qiao(xi, N2);
    X3 = x1.'*x2;
    [avgs(3, i), stds(3, i)] = comp_std(X3);
    X3 = abs(X3(:)).^2;
    paprs3(i) = max(X3)/mean(X3);
end

close all
figure
paprs = pow2db(paprs);
paprs3 = pow2db(paprs3);
plot(paprs, stds(1, :), '--r', LineWidth=1.5, Marker='o')
hold on
plot(paprs, stds(2, :), '--g', LineWidth=1.5, Marker='square')
plot(paprs3, stds(3, :), '--b', LineWidth=1.5, Marker='*')
legend('Opt wrt circle', 'Opt wrt square', 'Method in [12]', 'FontSize', 12)
xlabel('PAPR / dB')
ylabel('STD')
gca.FontSize = 12;

%% 
function [avg, STD] = comp_std(X)
[N1, N2] = size(X);
mode = 'uniform'; % the mode of choosing discrete angles
beamSampleHorizonNum = 360;
beamSampleVerticalNum = 180;
broadbeampattern = zeros(beamSampleHorizonNum, beamSampleVerticalNum);
[beamThetaVec, beamPhiVec] = gen_angle_vec(beamSampleHorizonNum, beamSampleVerticalNum, mode);
X = X/norm(X, 'fro');
for i = 1:beamSampleHorizonNum
    for j = 1:beamSampleVerticalNum
        u = -sin(beamPhiVec(j));
        v = -cos(beamThetaVec(i))*cos(beamPhiVec(j));
        F  = exp(1j*pi*(u*(0:N1-1).' + (v*(0:N2-1))));
        broadbeampattern(i, j) = abs(F(:).'*X(:))^2;
    end
end
avg = mean(broadbeampattern(:));
STD = std(broadbeampattern(:));
end




