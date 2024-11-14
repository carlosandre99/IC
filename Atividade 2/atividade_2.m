load('IN_OUT_PA.mat'); 

P = 5; 
M = 1; 

N = length(in);
XX = zeros(N, (P+1) * (M+1));

for n = 1:N
    col = 1;
    for p = 0:P
        for m = 0:M
            if n-m > 0
                XX(n, col) = in(n-m)^p;
            else
                XX(n, col) = 0;
            end
            col = col + 1;
        end
    end
end

H = XX \ out;
out_est = XX * H;

figure;
plot(out, 'b', 'DisplayName', 'Saída medida');
hold on;
plot(out_est, 'r', 'DisplayName', 'Saída estimada');
xlabel('Amostras');
ylabel('Saída');
legend;