clc; clear all; close all;
load('in_out_SBRT2_direto.mat');

% Modelagem inversa
% temp = in_extraction;
% in_extraction = out_extraction;
% out_extraction = temp;

fc = 900e6;
fc_mod = 900e6;
fs1 = 61.44e6;
fs2 = 8.97024e9;
alfa = 0;
phi = 0;
N_fft = 2048;

% Modulador/Demodulador ideal
fc_demod = 900e6;
[out_extraction_mod_ideal, out_validation_mod_ideal] = modulacao_demodulacao(out_extraction, out_validation, fc_mod, fc_demod, fs1, fs2, alfa, phi);

% Modulador/Demodulador não ideal com fc_demod = 910 MHz
fc_demod = 900e6;
omega1 = 2*pi*fc_mod;
omega2 = 2*pi*fc_demod;
t = (0:length(in_extraction)-1)'/fs1;

in_nao_ideal = in_extraction .* exp(1j*(omega1 - omega2)*t);
out_nao_ideal = out_extraction .* exp(1j*(omega1 - omega2)*t);

[out_extraction_mod_nao_ideal, out_validation_mod_nao_ideal] = modulacao_demodulacao(out_extraction, out_validation, fc_mod, fc_demod, fs1, fs2, alfa, phi);

% Caso Não Ideal: fc_demod = 910 MHz
signals_nao_ideal = mod_demod(out_extraction, fc_mod, fc_demod, fs1, fs2, alfa, phi);

% Plotting para o Caso Não Ideal
figure;
%sgtitle(sprintf('fc\\_demod = %.0f MHz', fc_demod / 1e6));

% Espectro do sinal original em 61.44 MHz
subplot(3, 2, 1);
espectro_frequencia_2(in_nao_ideal, fs1, '', 'b');
title('(a) Sinal de entrada a 61,44 MHz');

% Espectro do sinal reamostrado
subplot(3, 2, 2);
plot_spectrum(signals_nao_ideal.yi_reamostrado, fs2, false, 'b');
title('(b) Reamostragem 9 GHz');

% Espectro do sinal modulado (s_n)
subplot(3, 2, 3);
plot_spectrum(signals_nao_ideal.s_n, fs2, false, 'b');
title('(c) Sinal modulado s(n)');

% Espectro do sinal demodulado (antes da filtragem)
subplot(3, 2, 4);
plot_spectrum(signals_nao_ideal.yi_demod, fs2, false, 'b');
title('(d) Sinal demodulado (sem filtragem)');

% Espectro de saída RF (sinal filtrado)
subplot(3, 2, 5);
plot_spectrum(signals_nao_ideal.yi_filtrado, fs2, false, 'b');
title('(e) Sinal demodulado (filtrado)');

% Espectro final reamostrado para 61.44 MHz
subplot(3, 2, 6);
espectro_frequencia_2(out_extraction_mod_nao_ideal, fs1, '', 'b');
title('(f) Sinal de saída reamostado em 61,44 MHz');

%-------------------------------------------------------------------------%

% Plotting no domínio do tempo para as primeiras 200 amostras
figure;
%sgtitle(sprintf('fc\\_demod = %.0f MHz', fc_demod / 1e6));

% Sinal in_nao_ideal no tempo
subplot(3, 2, 1);
plot(t(1:200), real(in_nao_ideal(1:200)), 'b');
hold on;
plot(t(1:200), imag(in_nao_ideal(1:200)), 'r');
title('(a) Sinal de entrada a 61,44 MHz');
xlabel('Tempo (s)');
ylabel('Amplitude');
legend('Real', 'Imag');
%grid on;

% Sinal reamostrado yi_reamostrado no tempo 
t_resample = (0:length(signals_nao_ideal.yi_reamostrado)-1)'/fs2;
subplot(3, 2, 2);
plot(t_resample(1:min(3221, length(t_resample))), signals_nao_ideal.yi_reamostrado(1:min(3221, length(signals_nao_ideal.yi_reamostrado))), 'b');
title('(b) Reamostragem 9 GHz');
xlabel('Tempo (s)');
ylabel('Amplitude');
%grid on;

% Sinal modulado s_n no tempo 
t_s_n = (0:length(signals_nao_ideal.s_n)-1)'/fs2;
subplot(3, 2, 3);
plot(t_s_n(1:min(200, length(t_s_n))), signals_nao_ideal.s_n(1:min(200, length(signals_nao_ideal.s_n))), 'b');
title('(c) Sinal modulado s(n)');
xlabel('Tempo (s)');
ylabel('Amplitude');
%grid on;

% Sinal demodulado (sem filtragem) yi_demod no tempo
t_demod = (0:length(signals_nao_ideal.yi_demod)-1)'/fs2;
subplot(3, 2, 4);
plot(t_demod(1:min(200, length(t_demod))), signals_nao_ideal.yi_demod(1:min(200, length(signals_nao_ideal.yi_demod))), 'b');
title('(d) Sinal demodulado (sem filtragem)');
xlabel('Tempo (s)');
ylabel('Amplitude');
%grid on;

% Sinal demodulado filtrado no tempo
t_filtrado = (0:length(signals_nao_ideal.yi_filtrado)-1)'/fs2;
subplot(3, 2, 5);
plot(t_filtrado(1:min(3221, length(t_filtrado))), signals_nao_ideal.yi_filtrado(1:min(3221, length(signals_nao_ideal.yi_filtrado))), 'b');
title('(e) Sinal demodulado (filtrado)');
xlabel('Tempo (s)');
ylabel('Amplitude');
%grid on;

% Sinal final reamostrado para 61.44 MHz no tempo 
t_out = (0:length(out_extraction_mod_nao_ideal)-1)'/fs1;
subplot(3, 2, 6);
plot(t_out(1:min(200, length(t_out))), real(out_extraction_mod_nao_ideal(1:min(200, length(out_extraction_mod_nao_ideal)))), 'b');
hold on;
plot(t_out(1:min(200, length(t_out))), imag(out_extraction_mod_nao_ideal(1:min(200, length(out_extraction_mod_nao_ideal)))), 'r');
title('(f) Sinal de saída reamostado em 61,44 MHz');
xlabel('Tempo (s)');
ylabel('Amplitude');
legend('Real', 'Imag');
%grid on;

function plot_spectrum(signal, fs, use_log_scale, cor)
    N = length(signal);
    f = (-N/2:N/2-1)*(fs/N) / 1e9;  % Eixo de frequência em GHz
    spectrum = fftshift(abs(fft(signal))/N);  % Normalização 

    % Calcular a densidade de potência
    psd = (spectrum.^2) / (fs * N);  % Densidade de potência (W/Hz)
    
    % Converter para dBm/Hz
    psd_dBm = 10 * log10(psd * 1000);  % Convertendo W para mW e, em seguida, para dBm

    if use_log_scale
        semilogy(f, psd_dBm, cor); 
    else
        plot(f, psd_dBm, cor);
    end
    
    xlabel('Frequência (GHz)');  % Rótulo atualizado para GHz
    ylabel('PSD (dBm/Hz)');
  %  grid on;
end


function espectro_frequencia_2(sinal, fs, legenda, estilo)
    [Pxx, f] = pwelch(sinal, hanning(round(length(sinal)/8)), [], 4096, fs, 'twosided');
    Pxx_dBm = 10 * log10(Pxx * 1000);  

    f = f / 1e6;

    plot(f - fs/(2*1e6), fftshift(Pxx_dBm), estilo, 'DisplayName', legenda, 'LineWidth', 1.5); 
    xlabel('Frequência (MHz)');
    ylabel('PSD (dBm/Hz)');
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

function signals = mod_demod(out_extraction, fc_mod, fc_demod, fs1, fs2, alfa, phi)
    signals = struct();
    yi = real(out_extraction);
    yq = imag(out_extraction);
    
    % Reamostragem
    signals.yi_reamostrado = resample(yi, fs2, fs1);
    signals.yq_reamostrado = resample(yq, fs2, fs1);
    
    % Modulação
    t = (0:length(signals.yi_reamostrado)-1)'/fs2;
    signals.s_n = signals.yi_reamostrado .* cos(2*pi*fc_mod*t) - (1+alfa) * signals.yq_reamostrado .* sin(2*pi*fc_mod*t + phi);
    
    % Demodulação
    signals.yi_demod = 2 * signals.s_n .* cos(2*pi*fc_demod*t);
    signals.yq_demod = -2 * signals.s_n .* sin(2*pi*fc_demod*t);
    
    % Filtro passa-baixas
    N = min(length(signals.yi_demod), 2^nextpow2(length(signals.yi_demod)));
    f_filt = (-N/2:N/2-1)*(fs2/N);
    LPF = (abs(f_filt) < fc_mod); 
    
    yi_demod_freq = fftshift(fft(signals.yi_demod));
    yq_demod_freq = fftshift(fft(signals.yq_demod));
    
    signals.yi_filtrado = real(ifft(ifftshift(yi_demod_freq .* LPF')));
    signals.yq_filtrado = real(ifft(ifftshift(yq_demod_freq .* LPF')));
    
    % Reamostragem de volta para 61.44 MHz
    signals.yi_final = resample(signals.yi_filtrado, fs1, fs2);
    signals.yq_final = resample(signals.yq_filtrado, fs1, fs2);
end