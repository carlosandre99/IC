clc, clear all, close all;
load('in_out_SBRT2_direto.mat');

P = 3;
M = 1;
N = length(in_extraction);
XX = matriz_regressao(in_extraction, P, M, N);

% Modulador/Demodulador ideal
fc = 900e6;
fs1 = 61.44e6; 
fs2 = 8.97024e9;
alfa = 0;
phi = 0;
[out_extraction_mod_ideal, out_validation_mod_ideal] = modulacao_demodulacao(out_extraction, out_validation, fc, fc, fs1, fs2, alfa, phi);

H_ideal = XX \ out_extraction_mod_ideal;
N_val = length(in_validation);
XX_val = matriz_regressao(in_validation, P, M, N_val);
saida_estimada_ideal = XX_val * H_ideal;

nmse_ideal = 10 * log10(sum(abs(out_validation_mod_ideal - saida_estimada_ideal).^2) / sum(abs(out_validation_mod_ideal).^2));
disp(['NMSE com mod/demod ideal: ', num2str(nmse_ideal), ' dB']);

% Modulador/Demodulador não ideal com fc_demod = 910 MHz
fc_mod = 900e6;
fc_demod = 910e6;
[out_extraction_mod_naoideal, out_validation_mod_naoideal] = modulacao_demodulacao(out_extraction, out_validation, fc_mod, fc_demod, fs1, fs2, alfa, phi);

H_naoideal = [XX ones(N, 1)] \ out_extraction_mod_naoideal;
saida_estimada_naoideal = [XX_val ones(N_val, 1)] * H_naoideal;

nmse_naoideal_910 = 10 * log10(sum(abs(out_validation_mod_naoideal - saida_estimada_naoideal).^2) / sum(abs(out_validation_mod_naoideal).^2));
disp(['NMSE com mod/demod não ideal (fc_demod = 910 MHz): ', num2str(nmse_naoideal_910), ' dB']);

erro_freq_naoideal = abs(out_validation_mod_naoideal - saida_estimada_naoideal);
N_fft = 2048;
fs = 61.44e6; 
espectro_frequencia(erro_freq_naoideal, fs, N_fft, 'Erro (fc\_demod = 910 MHz)');

% Modulador/Demodulador não ideal com fc_demod = 901 MHz
fc_demod = 901e6;
[out_extraction_mod_naoideal, out_validation_mod_naoideal] = modulacao_demodulacao(out_extraction, out_validation, fc_mod, fc_demod, fs1, fs2, alfa, phi);

H_naoideal = [XX ones(N, 1)] \ out_extraction_mod_naoideal;
saida_estimada_naoideal = [XX_val ones(N_val, 1)] * H_naoideal;

nmse_naoideal_901 = 10 * log10(sum(abs(out_validation_mod_naoideal - saida_estimada_naoideal).^2) / sum(abs(out_validation_mod_naoideal).^2));
disp(['NMSE com mod/demod não ideal (fc_demod = 901 MHz): ', num2str(nmse_naoideal_901), ' dB']);

espectro_frequencia(in_validation, fs, N_fft, 'Espectro do sinal de entrada (in\_validation)');
espectro_frequencia(out_validation, fs, N_fft, 'Espectro do sinal de saída (out\_validation)');
espectro_frequencia(out_validation_mod_naoideal, fs, N_fft, 'Espectro da saída após demodulador não ideal');
espectro_frequencia(saida_estimada_ideal, fs, N_fft, 'Espectro da saída estimada pelo MP-1 (ideal)');
espectro_frequencia(saida_estimada_naoideal, fs, N_fft, 'Espectro da saída estimada pelo MP-2 (não ideal)');

erro_freq_naoideal = abs(out_validation_mod_naoideal - saida_estimada_naoideal);
espectro_frequencia(erro_freq_naoideal, fs, N_fft, 'Erro (fc\_demod = 901 MHz)');

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

function espectro_frequencia(sinal, fs, N_fft, titulo_texto)
    janela = hanning(length(sinal));
    sinal_janelado = sinal .* janela;
    espectro_freq = fftshift(fft(sinal_janelado, N_fft));
    f = (-N_fft/2:N_fft/2-1)*(fs/N_fft) / 1e6;
    figure;
    semilogy(f, abs(espectro_freq));
    title(titulo_texto);
    xlabel('Frequência (MHz)');
    ylabel('PSD');
end