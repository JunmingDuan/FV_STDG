function plot_ex2(K, n);

format long;
hold on;

%numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'.dat']);
%x1 = numer1(:,1); y1 = numer1(:,2);
%plot(x1, y1, 'o');

for K = 1:4;
  numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'.dat']);
  x1 = numer1(:,1); y1 = numer1(:,2);
  plot(x1, y1, '-');
end
legend('K=1', 'K=2', 'K=3', 'K=4');

