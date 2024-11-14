clc; clear all; close all;
load('in_out_SBRT2_direto.mat');

P = 3;  
M = 1;  
N = length(in_extraction);

in_extraction_abs = abs(in_extraction);
out_extraction_abs = abs(out_extraction);
in_validation_abs = abs(in_validation);
out_validation_abs = abs(out_validation);

XX = matriz_regressao(in_extraction_abs, P, M, N);

fc = 900e6;  
fs1 = 61.44e6;  
fs2 = 8.97024e9;  
alfa = 0;  
phi = 0;  

% Modulador/Demodulador não ideal com fc_demod = 910 MHz
fc_mod = 900e6;
fc_demod = 910e6;

[out_extraction_mod_naoideal, out_validation_mod_naoideal] = modulacao_demodulacao(out_extraction, out_validation, fc_mod, fc_demod, fs1, fs2, alfa, phi);

out_extraction_mod_naoideal_abs = abs(out_extraction_mod_naoideal);
out_validation_mod_naoideal_abs = abs(out_validation_mod_naoideal);

% Usando os dados de extração 
H_naoideal = [XX ones(N, 1)] \ out_extraction_mod_naoideal_abs;
saida_estimada_naoideal = [XX ones(N, 1)] * H_naoideal;

nmse_naoideal_extracao = 10 * log10(sum(abs(out_extraction_mod_naoideal_abs - saida_estimada_naoideal).^2) / sum(abs(out_extraction_mod_naoideal_abs).^2));
disp(['NMSE com mod/demod não ideal (dados de extração): ', num2str(nmse_naoideal_extracao), ' dB']);

abs_in_out(out_extraction_mod_naoideal_abs, in_extraction_abs, saida_estimada_naoideal, 'abs(out) x abs(in) - Não ideal (dados de extração)');
abs_out(out_extraction_mod_naoideal_abs, saida_estimada_naoideal, 'Não ideal (dados de extração)');
abs_erro(out_extraction_mod_naoideal_abs, saida_estimada_naoideal, 'Não ideal (dados de extração)');

% Usando dados de validação 
N_val = length(in_validation_abs);
XX_val = matriz_regressao(in_validation_abs, P, M, N_val);
saida_estimada_naoideal_val = [XX_val ones(N_val, 1)] * H_naoideal;

nmse_naoideal_validacao = 10 * log10(sum(abs(out_validation_mod_naoideal_abs - saida_estimada_naoideal_val).^2) / sum(abs(out_validation_mod_naoideal_abs).^2));
disp(['NMSE com mod/demod não ideal (dados de validação): ', num2str(nmse_naoideal_validacao), ' dB']);

disp('Coeficientes:');
disp(H_naoideal);

abs_in_out(out_validation_mod_naoideal_abs, in_validation_abs, saida_estimada_naoideal_val, 'abs(out) x abs(in) - Não ideal (dados de validação)');
abs_out(out_validation_mod_naoideal_abs, saida_estimada_naoideal_val, 'Não ideal (dados de validação)');
abs_erro(out_validation_mod_naoideal_abs, saida_estimada_naoideal_val, 'Não ideal (dados de validação)');

function abs_in_out(out_desejado, in_desejado, out_estimado, tipo)
    figure;
    plot(abs(in_desejado), abs(out_desejado), 'b.', 'DisplayName', 'Desejado');
    hold on;
    plot(abs(in_desejado), abs(out_estimado), 'r.', 'DisplayName', 'Estimado');
    title(tipo);
    xlabel('abs(in)');
    ylabel('abs(out)');
    legend;
end

function abs_out(out_validation, saida_estimada, tipo)
    t = (0:length(out_validation)-1);
    figure;
    plot(t, abs(out_validation), 'b', 'DisplayName', 'Desejado');
    hold on;
    plot(t, abs(saida_estimada), 'r', 'DisplayName', 'Estimado');
    title(['abs(out) - ', tipo]);
    xlabel('Amostra');
    ylabel('Magnitude');
    legend;
end

function abs_erro(out_validation, saida_estimada, tipo)
    t = (0:length(out_validation)-1);
    erro = abs(out_validation - saida_estimada);
    figure;
    plot(t, erro);
    title(['abs(erro) - ', tipo]);
    xlabel('Amostra');
    ylabel('Erro');
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
        LPF = (abs(f) < fc_mod);
        yi_demod_freq = fftshift(fft(yi_demod));
        yq_demod_freq = fftshift(fft(yq_demod));
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