load('in_out_SBRT2_direto.mat');

% Separar as componentes real e imaginária
yi = real(in_extraction);
yq = imag(in_extraction);

fc = 900e6; 
fs1 = 61.44e6; 
fs2 = 8.97024e9;  
alpha = 0; 
phi = 0; 

% Reamostragem do sinal
yi_resampled = resample(yi, fs2, fs1);
yq_resampled = resample(yq, fs2, fs1);

% Modulação de acordo com a equação (8)
t = (0:length(yi_resampled)-1)'/fs2;
s_n = yi_resampled .* cos(2*pi*fc*t) - (1+alpha) * yq_resampled .* sin(2*pi*fc*t + phi);

% Demodulação de acordo com a equação (9)
yi_demod = 2 * s_n .* cos(2*pi*fc*t);
yq_demod = -2 * s_n .* sin(2*pi*fc*t);

% Quantidade de pontos seja uma potência de 2
N = min(length(yi_demod), 2^nextpow2(length(yi_demod))); 
yi_demod = yi_demod(1:N);  
yq_demod = yq_demod(1:N);

f = (-N/2:N/2-1)*(fs2/N);
LPF = (abs(f) < fc); 

yi_demod_freq = fftshift(fft(yi_demod));
yq_demod_freq = fftshift(fft(yq_demod));

yi_filtered = real(ifft(ifftshift(yi_demod_freq .* LPF')));
yq_filtered = real(ifft(ifftshift(yq_demod_freq .* LPF')));

% Reamostragem de volta para 61.44 MHz
yi_final = resample(yi_filtered, fs1, fs2);
yq_final = resample(yq_filtered, fs1, fs2);

% Ajuste de tamanho dos sinais para cálculo e plotagem
min_length = min(length(yi), length(yi_final));
yi = yi(1:min_length);
yq = yq(1:min_length);
yi_final = yi_final(1:min_length);
yq_final = yq_final(1:min_length);

% Cálculo do NMSE 
yref = yi + 1i*yq;
ytest = yi_final + 1i*yq_final;
e_total = yref - ytest;
NMSE_total = 10*log10(sum(abs(e_total).^2) / sum(abs(yref).^2));
disp(['NMSE: ', num2str(NMSE_total), ' dB']);

figure;
subplot(2, 1, 1);
plot(yi, 'b', 'LineWidth', 2);
hold on;
plot(yi_final, 'r--', 'LineWidth', 1.5);
xlabel('Amostras');
ylabel('Parte real');
legend('Entrada', 'Saída');

subplot(2, 1, 2);
plot(yq, 'b', 'LineWidth', 2);
hold on;
plot(yq_final, 'r--', 'LineWidth', 1.5);
xlabel('Amostras');
ylabel('Parte Imaginária');
legend('Entrada', 'Saída');

N_fft = 2048;

figure;
subplot(2,1,1);
espectro_frequencia(yi, fs1, N_fft, 'b', 'Entrada');
hold on;
espectro_frequencia(yi_final, fs1, N_fft, 'r', 'Saída');
hold off;
title('Parte real');
legend;

subplot(2,1,2);
espectro_frequencia(yq, fs1, N_fft, 'b', 'Entrada');
hold on;
espectro_frequencia(yq_final, fs1, N_fft, 'r', 'Saída');
hold off;
title('Parte Imaginária');
legend;

function espectro_frequencia(sinal, fs, N_fft, cor, legenda)
    % Janela Hanning
    janela = hanning(length(sinal));        
    sinal_janelado = sinal .* janela;        
    
    % Calcula FFT 
    espectro_freq = fftshift(fft(sinal_janelado, N_fft));  
    
    % Calcula PSD em Watts/Hz
    Pxx = (abs(espectro_freq).^2) / (fs * sum(janela.^2)); 
    
    % Converte PSD para miliwatts (mW)
    Pxx_mW = Pxx * 1000;  
    
    % Converte para dBm/Hz
    Pxx_dBm = 10 * log10(Pxx_mW);
    
    % Vetor de frequências em MHz
    f = (-N_fft/2:N_fft/2-1)*(fs/N_fft) / 1e6;  
    
    plot(f, Pxx_dBm, cor, 'DisplayName', legenda);        
    xlabel('Frequência (MHz)');                                 
    ylabel('PSD (dBm/Hz)');                            
end
