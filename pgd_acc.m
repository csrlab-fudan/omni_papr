function [W, losses] = pgd_acc(M, N, k, papr, lp, epsilon, init, domain)
if nargin<9
    domain = 'circle';
end

M1 = M*k; N1 = N*k;
if isequal(domain, 'square')
    idx = ones(M1, N1)==1;
else
    alpha = 1.025;
    [indexCol, indexRow] = meshgrid(-1:2/M1:1-2/M1, -1:2/N1:1-2/N1);
    idx = (indexRow.^2 + indexCol.^2) <= alpha^2;
    idx = ifftshift(idx);
end

if isempty(init)
    w1 = exp(1j*pi/N*(0:M-1).*(1:M));
    w2 = exp(1j*pi/N*(0:M-1).*(1:M));
    W = w1.'*w2;
else
    W = init;
end


i = 1;
t = 1;
y = W;
s = 1;
N = 10;
iter_num = 5e4;
losses = zeros(1, iter_num);
fOld = cost_function(W, idx, lp);
deltas = zeros(1, N);
criterion = inf;
while criterion > epsilon && i<=iter_num
    WOld = W;
    [s, W] = armijo_search(y, idx, s, lp, lb, papr);
    tOld = t;
    t = 1/2*(1+sqrt(1+4*t^2));
    y = W+(tOld-1)/t*(W-WOld);
    f = cost_function(W, idx, lp);
    delta = abs(f-fOld)/fOld;
    deltas(mod(i-1, N)+1) = delta;
    criterion = mean(deltas(deltas~=0));
    fOld = f;
    fprintf('Iteration: %d; loss: %.6f; delta: %.5E\n', i, f, delta);
    losses(i) = f;
    i = i+1;
end
close all
figure
semilogy(losses(losses~=0))
end
%% cost function and gradient
function [f, g] = cost_function(W, idx, lp)
[M, N] = size(W);
[M1, N1] = size(idx);
% objective function
F = 1/sqrt(M*N) * fft2(W, M1, N1);
res = zeros(M1, N1);
if isequal(lp, 'ln')
    res(idx) =  2*log(abs(F(idx)));
else  
    res(idx) =  (abs(F(idx)).^lp)-1;
end
f = norm(res,"fro")^2;
% gradient
if nargout > 1
    if isequal(lp, 'ln')
        X = res.*F./(abs(F).^2);
    else
        X =  res.*F.*(lp/2*abs(F).^(lp-2));
    end
    g= 2*(M1*N1)/sqrt(M*N) * ifft2(X, M1, N1);
    g = g(1:M, 1:N);
end
end
%% armijo
function [s, x] = armijo_search(y, idx, s, lp, papr)
beta = 1/2;
[fOld, g] = cost_function(y, idx, lp);
while  true  % to find s
    x = project(y - s*2*g, papr);
    fNew = cost_function(x, idx, lp);
    delta = reshape((x-y), 1, []);
    if fNew-fOld > real(delta*2*conj(g(:)))+1/(2*s)*norm(delta)^2
        s = beta*s;
        continue;
    end
    break;
end
end
%% projection
function Wp = project(W, p)
[M, N] = size(W);
r = sort(abs(W(:)), 'descend');
j = 0;
while true
    eta = sqrt(sum(r(j+1:end).^2)/(M*N-j*p));
    if r(j+1) > eta*sqrt(p)
        j = j+1;
    else
        break
    end
end
rho = abs(W);
rhox = sqrt(p)*ones(M, N);
index = rho/eta <= sqrt(p);
rhox(index) = rho(index)/eta;
Wp = rhox .* exp(1j*angle(W));
end