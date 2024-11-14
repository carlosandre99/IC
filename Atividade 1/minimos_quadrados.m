Y = [0.2; 0.3; 0.45; 0.7; 0.8];
X = [0; 0.1; 0.2; 0.3; 0.4];
XX = [X, ones(size(X))];

COEF = XX \ Y;
a = COEF(1);
b = COEF(2);

fprintf('Coeficientes ajustáveis: a = %.4f, b = %.4f\n', a, b);

Y_est = XX * COEF;

figure;
scatter(X, Y, 'b', 'filled');
hold on;
plot(X, Y_est, 'r');
xlabel('xi');
ylabel('yi');
legend('Pontos', 'Reta ajustada');
