function plot_ex1(K, n);

format long;
hold on;
mu = 1;
t = 0.5;

%f = @(x) sin(2*pi*(x - mu*t))*exp(-t);
f = @(x) sin(2*pi*(x - mu*t));
%f = @(x) sin(2*pi*(x))*exp(-t);
%f = @(x) exp(-1e-1*(t*ones(size(x))));
x0 = linspace(0, 1, 1e2);
plot(x0, f(x0), '-r');

numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'.dat']);
x1 = numer1(:,1); y1 = numer1(:,2);
plot(x1, y1, 'o');

m = 6;
err = zeros(m, 3);
for i = 1:m;
  n = 10*2^(i-1);
  %n = 2^(i+1);
  numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'.dat']);
  x1 = numer1(:,1); y1 = numer1(:,2);
  err(i, 1) = sqrt(sum((y1-f(x1)).^2)/n);
  err(i, 2) = sum(abs(y1-f(x1)))/n;
  err(i, 3) = max(abs(y1-f(x1)));
end
order = zeros(m-1, 3);
for i = 1:m-1;
  order(i,:) = log2(err(i,:)./err(i+1,:));
end
err
order

%for K = 1:4;
  %numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'.dat']);
  %x1 = numer1(:,1); y1 = numer1(:,2);
  %plot(x1, y1, '-');
%end
%legend('K=1', 'K=2', 'K=3', 'K=4');

