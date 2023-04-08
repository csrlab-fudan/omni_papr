function [x, papr] = gen_qiao(xi, M)
theta = asin(2/(2*M-1)*((1:2*M-1)-M));
epsilon = xi*sin(theta);
W = exp(1j*pi*(1-M:M-1).'*sin(theta));
P = inv(W);
r = 1+epsilon;
a = r*P(:, 1:M);
a = [conj(a(1:end-1)), flip(a)];
solutions = roots(a);
[~, idx] = sort(abs(solutions));
solutions = solutions(idx);
[x, papr] = chose_pair(solutions);
x = x*sqrt(M);
% figure
% semilogy(abs(xcorr(x)));
% hold on
% semilogy(abs(a));
% legend("corr", "a")

papr = pow2db(papr);
end

function [y, p] = chose_pair(solutions)
N = length(solutions)/2+1;
ds = de2bi(0:2^(N-1)-1);
p = inf;
for i = 1:size(ds, 1)
    d = ds(i, :);
    a = zeros(1, N-1);
    for j = 1:N-1
        if d(j)==0
            a(j) = solutions(j);
        else
            a(j) = 1/conj(solutions(j));
        end
    end
    x = poly(a);
    papr = N*max(abs(x).^2)/norm(x)^2;
    if papr < p
        y = x/norm(x);
        p = papr;
    end
end
end
