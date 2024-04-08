Ns = 10.^(1:6);
k = 2;
fxs = cell(size(Ns));
vars = zeros(size(Ns));
for ii = 1:length(Ns)
    ii
    N = Ns(ii);
    idx = -1:2/(k*N):1-2/(k*N);
    x = exp(1j*pi/N*(0:N-1).*(1:N)).';
    x = x/norm(x);
    fx = fftshift(fft(x, k*N));
    fx = abs(fx).^2;
    avgii = 0;
    varii = 0;
    L = 0;
    for i = 1:k*N
        if mod(i, k*N/100)==1
            (i-1)/(k*N/100)
        end
        inCircle = idx(i)^2+idx.^2<=1;
        temp = fx(i)*fx(inCircle);
        L = L + length(temp);
        varii = varii + sum((temp-1).^2);
    end
    vars(ii) = 1/L*varii;
end
close all
figure
loglog(Ns, vars, '--r', 'LineWidth', 1.5, Marker='*')
hold on
anal = 2/pi*ones(size(Ns))./(Ns.^0.5);
anal = 2*anal + anal.*anal;
loglog(Ns, anal, '--g', 'LineWidth', 1.5, Marker='o')
xlabel('N', 'FontSize',12)
ylabel('NMSE', 'FontSize',12)
legend('Simulation', 'Analysis', 'FontSize', 12)





