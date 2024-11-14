clc; clear all; close all;
load('in_out_SBRT2_direto.mat');

P = 4; 
M = 4;
N = length(in_extraction);
N_val = length(in_validation);

fc_mod = 900e6;
fc_demod = 915e6;
fs1 = 61.44e6;

omega1 = 2*pi*fc_mod;
omega2 = 2*pi*fc_demod;

t_extraction = (0:N-1)'/fs1;
t_validation = (0:N_val-1)'/fs1;

in_nao_ideal = in_extraction .* exp(1j*(omega1 - omega2)*t_extraction);
out_nao_ideal = out_extraction .* exp(1j*(omega1 - omega2)*t_extraction);

% Valores máximos das partes reais e imaginárias 
max_real_in_extraction = max(abs(real(in_extraction)));
max_real_out_extraction = max(abs(real(out_extraction)));
max_imag_in_extraction = max(abs(imag(in_extraction)));
max_imag_out_extraction = max(abs(imag(out_extraction)));
max_real_in_validation = max(abs(real(in_validation)));
max_real_out_validation = max(abs(real(out_validation)));
max_imag_in_validation = max(abs(imag(in_validation)));
max_imag_out_validation = max(abs(imag(out_validation)));

% Maior dos máximos
maior_maximo = max([max_real_in_extraction, max_real_out_extraction, max_imag_in_extraction, max_imag_out_extraction, ...
                    max_real_in_validation, max_real_out_validation, max_imag_in_validation, max_imag_out_validation]);

% Normalização 
in_extraction = in_extraction / maior_maximo;
out_extraction = out_extraction / maior_maximo;
in_validation = in_validation / maior_maximo;
out_validation = out_validation / maior_maximo;

% Parte 1 - Treinamento do modelo MP do PA
XX = matriz_regressao(in_nao_ideal, P, M, N);
H = XX \ out_nao_ideal;

XX_val = matriz_regressao(in_validation, P, M, N_val);
out_estimado_PA = XX_val * H;

out_estimado_PA = out_estimado_PA .* exp(1j*(omega1 - omega2)*t_validation);

% Parte 2 - Treinamento do modelo MP do DPD 
XX_DPD = matriz_regressao(out_nao_ideal, P, M, N);
H_DPD = XX_DPD \ in_nao_ideal;

ganho = abs(out_estimado_PA) ./ abs(in_validation); 
in_cascata = ganho .* in_validation .* exp(1j*(omega1 - omega2)*t_validation);

XX_val_DPD = matriz_regressao(in_cascata, P, M, N_val);
out_estimado_DPD = XX_val_DPD * H_DPD;

% Parte 3 - Cascata
XX_val_PA = matriz_regressao(out_estimado_DPD, P, M, N_val);
out_cascata = XX_val_PA * H;

nmse_cascata = 10 * log10(sum(abs(in_cascata - out_cascata).^2) / sum(abs(in_cascata).^2));
disp(['NMSE da cascata: ', num2str(nmse_cascata), ' dB']);

% AM-AM
figure;
hold on;
plot(abs(out_estimado_DPD), abs(out_estimado_PA), 'b.', 'DisplayName', 'Curva PA');
plot(abs(in_cascata), abs(out_estimado_DPD), 'r.', 'DisplayName', 'Curva DPD');
plot(abs(in_cascata), abs(out_cascata), 'g.', 'DisplayName', 'Curva cascata (DPD + PA)');
xlabel('|in|');
ylabel('|out|');
title(['AM-AM ', num2str(fc_demod/1e6), ' MHz']);  
legend('Location', 'Best');

% AM-PM
figure;
hold on;
plot(abs(in_cascata), angle(out_estimado_PA) - angle(in_cascata), 'b.', 'DisplayName', 'PA'); 
plot(abs(in_cascata), angle(out_estimado_DPD) - angle(in_cascata), 'r.', 'DisplayName', 'DPD');
plot(abs(in_cascata), angle(out_cascata) - angle(in_cascata), 'g.', 'DisplayName', 'Saída cascata (DPD + PA)');
xlabel('|in|');
ylabel('Diferença de fase');
title(['AM-PM ', num2str(fc_demod/1e6), ' MHz']);  
legend;

ganho = mean(abs(out_estimado_PA)) ./ mean(abs(in_validation)); 
in_cascata = ganho .* in_validation .* exp(1j*(omega1 - omega2)*t_validation);

XX_val_DPD = matriz_regressao(in_cascata, P, M, N_val);
out_estimado_DPD = XX_val_DPD * H_DPD;

XX_val_PA = matriz_regressao(out_estimado_DPD, P, M, N_val);
out_cascata = XX_val_PA * H;

% PSD
fs1 = 61.44e6;
figure;
hold on;
espectro_frequencia_2(in_cascata, fs1, 'Entrada cascata', 'b'); 
%espectro_frequencia_2(out_estimado_DPD, fs1, 'Saída DPD', 'r');
espectro_frequencia_2(out_estimado_PA, fs1, 'Saída PA (sem DPD)', 'g'); 
espectro_frequencia_2(out_cascata, fs1, 'Saída cascata (DPD + PA)', 'k');
ylabel('Densidade Espectral de Potência (dBm/Hz)');  
legend;
hold off;

% Tempo
figure;
subplot(2,1,1);
plot(real(in_cascata), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Entrada da cascata');
hold on;
plot(real(out_cascata), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Saída da cascata');
xlabel('Amostra');
ylabel('Parte Real');  
ylim([-1 1]);
legend('Location', 'Best');

subplot(2,1,2);
plot(imag(in_cascata), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Entrada da cascata');
hold on;
plot(imag(out_cascata), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Saída da cascata');
xlabel('Amostra');
ylabel('Parte Imaginária');
ylim([-1 1]);
legend('Location', 'Best');

function espectro_frequencia_2(sinal, fs, legenda, estilo)
    [Pxx, f] = pwelch(sinal, hanning(round(length(sinal)/8)), [], 4096, fs, 'twosided');
    Pxx_dBm = 10 * log10(Pxx * 1000);  

    f = f / 1e6;

    plot(f - fs/(2*1e6), fftshift(Pxx_dBm), estilo, 'DisplayName', legenda, 'LineWidth', 1.5); 
    xlabel('Frequência (MHz)');
    ylabel('Densidade Espectral de Potência (dBm/Hz)');
end

function XX = matriz_regressao(sinal, P, M, N)
    XX = zeros(N, P * (M+1));
    for n = 1:N
        coluna = 1;
        for p = 1:P
            for m = 0:M
                if n-m > 0
                    XX(n, coluna) = abs(sinal(n-m))^(2*(p-1)) * sinal(n-m);
                else
                    XX(n, coluna) = 0;
                end
                coluna = coluna + 1;
            end
        end
    end
end
