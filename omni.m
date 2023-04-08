clear;
antRow = 16;
antCol = 16;
papr = db2pow(3);
scaleFactor = 1;   % the circle scale factor
extendFactor = 1.025; % the circle extending factor to diminish the classification error
overSamplingRate = 10;
N1 = antRow; M1= overSamplingRate*N1;
N2 = antCol; M2 = overSamplingRate*N2;


beamSampleHorizonNum = 360;
beamSampleVerticalNum = 180;
beamSampleLength = beamSampleHorizonNum*beamSampleVerticalNum;
broadbeampattern = zeros(beamSampleHorizonNum, beamSampleVerticalNum);
mode = 'uniform'; % the mode of choosing discretized angles
[beamThetaVec, beamPhiVec] = gen_angle_vec(beamSampleHorizonNum, beamSampleVerticalNum, mode);


epsilon = 1e-8;
lp = 1; % lp norm
X1 = pgd_acc(N1, N2, overSamplingRate, papr, lp, epsilon, []);
X2 = pgd_acc(N1, N2, overSamplingRate, papr, lp, epsilon, [], 'square');

xi = 0.02;
[x1, ~] = gen_qiao(xi, N1);
[x2, ~] = gen_qiao(xi, N2);
X3 = x1.'*x2;
Xs = cat(3, X1, X2, X3);


%%
close all
figure
t = tiledlayout('flow');
linetype = {'-', '--', '-.'};
mark = {'o', 's'};
nexttile
for ii = 1:3
    X = Xs(:, :, ii);
    h = cdfplot(abs(X(:)).^2);
    h.LineStyle = linetype{ii};
    h.LineWidth = 1.5;
    hold on
end
xlabel('Transmitted power', 'FontSize',12);
ylabel('CDF', 'FontSize',12)
legend('Circle', 'Square', 'Method in [12]', 'FontSize', 12)
title('')

nexttile
for ii = 1:3
    X = Xs(:, :, ii);
    X = X/norm(X, 'fro');
    for i = 1:beamSampleHorizonNum
        for j = 1:beamSampleVerticalNum
            u = -sin(beamPhiVec(j));
            v = -cos(beamThetaVec(i))*cos(beamPhiVec(j));
            F  = exp(1j*pi*scaleFactor*(u*(0:N1-1).' + (v*(0:N2-1))));
            broadbeampattern(i, j) = abs(F(:).'*X(:))^2;
        end
    end
    h = cdfplot(broadbeampattern(:));
    h.LineStyle = linetype{ii};
    h.LineWidth = 1.5;
    hold on
end
title('')
xlabel('Received power', 'FontSize',12);
ylabel('CDF', 'FontSize',12)
legend('Circle', 'Square', 'Method in [9]', 'FontSize',12)
t.TileSpacing = 'compact';
t.Padding = 'compact';
    

