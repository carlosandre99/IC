load('in_out_SBRT2_direto.mat');

P = 3; 
M = 1;

N = length(in_extraction);

% Matriz de regressão XX para o modelo MP - Sinais de extração
XX = zeros(N, P * (M+1));
for n = 1:N
    col = 1;
    for p = 1:P
        for m = 0:M
            if n-m > 0
                XX(n, col) = abs(in_extraction(n-m))^(2*(p-1)) * in_extraction(n-m);
            else
                XX(n, col) = 0;
            end
            col = col + 1;
        end
    end
end

H = XX \ out_extraction;

% Saída estimada pelo modelo 

N_val = length(in_validation);

% Matriz de regressão XX para o modelo MP - Sinais de de validação
XX_val = zeros(N_val, P * (M+1));
for n = 1:N_val
    col = 1;
    for p = 1:P
        for m = 0:M
            if n-m > 0
                XX_val(n, col) = abs(in_validation(n-m))^(2*(p-1)) * in_validation(n-m);
            else
                XX_val(n, col) = 0;
            end
            col = col + 1;
        end
    end
end

out_est_val = XX_val * H;

nmse = 10 * log10(sum(abs(out_validation - out_est_val).^2) / sum(abs(out_validation).^2));
disp(['NMSE: ', num2str(nmse), ' dB']);

figure;
plot(abs(in_validation), abs(out_validation), 'b.', 'DisplayName', 'Medido');
hold on;
plot(abs(in_validation), abs(out_est_val), 'r.', 'DisplayName', 'Estimado');
xlabel('Amplitude de entrada');
ylabel('Amplitude de saída');
legend;
title('Curva AM-AM');

figure;
plot(abs(in_validation), angle(out_validation) - angle(in_validation), 'b.', 'DisplayName', 'Medido');
hold on;
plot(abs(in_validation), angle(out_est_val) - angle(in_validation), 'r.', 'DisplayName', 'Estimado');
xlabel('Amplitude de entrada');
ylabel('Diferença de fase');
legend;
title('Curva AM-PM');