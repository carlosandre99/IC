clc; clear all; close all;
load('in_out_SBRT2_direto.mat');

% Troca de in_extraction e out_extraction
% temp = in_extraction;
% in_extraction = out_extraction;
% out_extraction = temp;

P = 4; 
M = 4;
N = length(in_extraction);

fc = 900e6;
fc_mod = 900e6;
fs1 = 61.44e6;
fs2 = 8.97024e9;
alfa = 0;
phi = 0;

% Modulador/Demodulador ideal
fc_demod = 900e6;

[out_extraction_mod_ideal, out_validation_mod_ideal] = modulacao_demodulacao(out_extraction, out_validation, fc_mod, fc_demod, fs1, fs2, alfa, phi);

XX_ideal = matriz_regressao(in_extraction, P, M, N);
H_ideal = XX_ideal \ out_extraction;
saida_estimada_ideal = XX_ideal * H_ideal;

nmse_ideal = 10 * log10(sum(abs(out_extraction_mod_ideal - saida_estimada_ideal).^2) / sum(abs(out_extraction_mod_ideal).^2));
disp(['NMSE para o caso ideal (fc_demod = ', num2str(fc_demod/1e6), ' MHz): ', num2str(nmse_ideal), ' dB']);

% Modulador/Demodulador não ideal
fc_demod = 910e6;
omega1 = 2*pi*fc_mod;
omega2 = 2*pi*fc_demod;
t = (0:length(in_extraction)-1)'/fs1;

in_nao_ideal = in_extraction .* exp(1j*(omega1 - omega2)*t);
out_nao_ideal = out_extraction .* exp(1j*(omega1 - omega2)*t);

[out_extraction_mod_nao_ideal, out_validation_mod_nao_ideal] = modulacao_demodulacao(out_extraction, out_validation, fc_mod, fc_demod, fs1, fs2, alfa, phi);

XX_nao_ideal = matriz_regressao(in_nao_ideal, P, M, N);
H_nao_ideal = XX_nao_ideal \ out_nao_ideal;
saida_estimada_nao_ideal = XX_nao_ideal * H_nao_ideal;

nmse_nao_ideal = 10 * log10(sum(abs(out_extraction_mod_nao_ideal - saida_estimada_nao_ideal).^2) / sum(abs(out_extraction_mod_nao_ideal).^2));
disp(['NMSE para o caso não ideal (fc_demod = ', num2str(fc_demod/1e6), ' MHz): ', num2str(nmse_nao_ideal), ' dB']);

%-------------------------------------------------------------------------%
N_fft = 2048;

figure;
subplot(2,1,1);
espectro_frequencia(in_extraction, fs1, N_fft, 'Extração');
hold on;
espectro_frequencia(in_nao_ideal, fs1, N_fft, 'Não ideal');
hold off;
title('Sinais de entrada');
ylabel('|in|');  
legend;

subplot(2,1,2);
espectro_frequencia(out_extraction, fs1, N_fft, 'Extração');
hold on;
espectro_frequencia(out_extraction_mod_nao_ideal, fs1, N_fft, 'Não ideal');
espectro_frequencia(saida_estimada_nao_ideal, fs1, N_fft, 'Estimada');
hold off;
title('Sinais de saída');
ylabel('|out|'); 
legend;
%-------------------------------------------------------------------------%

% AM-PM
figure;
hold on;
plot(abs(in_nao_ideal), angle(out_extraction_mod_nao_ideal) - angle(in_nao_ideal), 'b.', 'DisplayName', 'Desejado');
plot(abs(in_nao_ideal), angle(saida_estimada_nao_ideal) - angle(in_nao_ideal), 'r.', 'DisplayName', 'Estimado');
xlabel('|in|');
ylabel('Diferença de fase');
title(['AM-PM (fc\_demod = ', num2str(fc_demod/1e6), ' MHz)']);  
legend;

% abs(out) x abs(in) - Saída MP
figure;
hold on;
plot(abs(in_nao_ideal), abs(out_extraction_mod_nao_ideal), 'bo', 'DisplayName', 'Curva PA');
plot(abs(in_nao_ideal), abs(saida_estimada_nao_ideal), 'ro', 'DisplayName', 'Curva estimada');
xlabel('|in|');
ylabel('|out|');
%title(['AM-AM (fc\_demod = ', num2str(fc_demod/1e6), ' MHz)']);  
legend('Location', 'Best');

% Primeiras 200 amostras
figure;
subplot(2,1,1);
hold on;
plot(t(1:200), real(out_extraction_mod_nao_ideal(1:200)), 'b-', 'DisplayName', 'Saída PA');
plot(t(1:200), real(saida_estimada_nao_ideal(1:200)), 'r*', 'DisplayName', 'Saída estimada');
hold off;
title('Parte real');
xlabel('Tempo (s)');
ylabel('Amplitude ');
legend('Location', 'Best');

subplot(2,1,2);
hold on;
plot(t(1:200), imag(out_extraction_mod_nao_ideal(1:200)), 'b-', 'DisplayName', 'Saída PA');
plot(t(1:200), imag(saida_estimada_nao_ideal(1:200)), 'r*', 'DisplayName', 'Saída estimada');
hold off;
title('Parte imaginária');
xlabel('Tempo (s)');
ylabel('Amplitude');
legend('Location', 'Best');

function espectro_frequencia(sinal, fs, N_fft, legenda)
    janela = hanning(length(sinal));        
    sinal_janelado = sinal .* janela;        
    espectro_freq = fftshift(fft(sinal_janelado, N_fft));  
    f = (-N_fft/2:N_fft/2-1)*(fs/N_fft) / 1e6;  
    semilogy(f, abs(espectro_freq), 'DisplayName', legenda);        
    xlabel('Frequência (MHz)');                                                               
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

function [out_extraction_mod, out_validation_mod] = modulacao_demodulacao(out_extraction, out_validation, fc_mod, fc_demod, fs1, fs2, alfa, phi)
    [out_extraction_mod, out_validation_mod] = deal(zeros(size(out_extraction)));
    for idx = 1:2
        if idx == 1
            y = out_extraction;
        else
            y = out_validation;
        end
        yi = real(y);
        yq = imag(y);
        
        % Reamostragem
        yi_reamostrado = resample(yi, fs2, fs1);
        yq_reamostrado = resample(yq, fs2, fs1);
        
        % Modulação
        t = (0:length(yi_reamostrado)-1)'/fs2;
        s_n = yi_reamostrado .* cos(2*pi*fc_mod*t) - (1+alfa) * yq_reamostrado .* sin(2*pi*fc_mod*t + phi);
        
        % Demodulação
        yi_demod = 2 * s_n .* cos(2*pi*fc_demod*t);
        yq_demod = -2 * s_n .* sin(2*pi*fc_demod*t);
        N = min(length(yi_demod), 2^nextpow2(length(yi_demod)));
        yi_demod = yi_demod(1:N);
        yq_demod = yq_demod(1:N);
        f = (-N/2:N/2-1)*(fs2/N);
        
        % Filtro passa-baixas
        LPF = (abs(f) < fc_mod); 
        yi_demod_freq = fftshift(fft(yi_demod));
        yq_demod_freq = fftshift(fft(yq_demod));
        
        % Aplicando o filtro
        yi_filtrado = real(ifft(ifftshift(yi_demod_freq .* LPF')));
        yq_filtrado = real(ifft(ifftshift(yq_demod_freq .* LPF')));
        
        % Reamostragem de volta para 61.44 MHz
        yi_final = resample(yi_filtrado, fs1, fs2);
        yq_final = resample(yq_filtrado, fs1, fs2);
        
        % Ajuste de tamanho dos sinais
        min_length = min(length(yi), length(yi_final));
        yi_final = yi_final(1:min_length);
        yq_final = yq_final(1:min_length);
        
        if idx == 1
            out_extraction_mod = yi_final + 1i*yq_final;
        else
            out_validation_mod = yi_final + 1i*yq_final;
        end
    end
end