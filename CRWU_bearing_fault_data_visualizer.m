% Ladataan data
load('CWRU Inner 1730 21.mat'); 

% Oletetaan, että muuttujan nimi on X234_FE_time (fan end data)
signal = X212_DE_time;

% Parametrit
Fs = 12000;          % Näytteenottotaajuus (Hz)
N = length(signal);  % Näytteiden määrä
t = (0:N-1)/Fs;       % Aikavektori

% Määritetään haluttu kesto 
duration_sec = 0.4;
samples_to_plot = round(duration_sec * Fs);

% Rajataan signaali ja aika
t_short = t(1:samples_to_plot);
signal_short = signal(1:samples_to_plot);

% Laskennallinen BPFO ja BPFI
N_b = 9;
D_b = 0.3126;
D_c = 1.537;
Angle = 0;
Fr = X212RPM/60;
BPFO = N_b/2*Fr*(1-D_b*cos(Angle)/D_c);
BPFI = N_b/2*Fr*(1+D_b*cos(Angle)/D_c);
BSF = D_c/(2*D_b)*Fr*(1-(D_b*cos(Angle)/D_c)^2);
FTF = 1/2*Fr*(1-(D_b*cos(Angle)/D_c));
F1_mod_positive = Fr + BPFI;
F1_mod_negative = Fr - BPFI;
F2_mod_positive = FTF + BPFI;


%% Piirretään aikatason signaali (rajattu)
figure;
subplot(3,1,1);       
plot(t_short, signal_short);
title('Aikatason signaali (Käyttöpää)');
xlabel('Aika (s)');
ylabel('Kiihtyvyys');
grid on;

%% FFT
window = hamming(N);
signal_windowed = signal .* window;
Y = fft(signal_windowed);
half_N = floor(N/2);
f = Fs*(0:half_N)/N;        
P2 = abs(Y/N);              
P1 = P2(1:half_N+1);        
P1(2:end-1) = 2*P1(2:end-1); 

% Rajataan taajuusalue
max_freq = 550;
f_range_idx = f <= max_freq;
f_plot = f(f_range_idx);
P1_plot = P1(f_range_idx);


% Etsitään piikit
[peaks, locs] = findpeaks(P1_plot, 'MinPeakHeight', 0.004, ... 
                          'MinPeakDistance', round(10 * length(P1) / Fs));
peak_freqs = f_plot(locs);
peak_amps = P1_plot(locs);

% Muutetaan amplitudi desibeleiksi
P1_dB = 20*log10(P1_plot);

% Piirretään FFT-spektri logaritmisella skaalaalla 
subplot(3,1,2);
plot(f_plot, P1_dB, 'LineWidth',1.5);
hold on;
ylim([-120, 0]); 
plot(peak_freqs, 20*log10(peak_amps), 'ro', 'MarkerFaceColor', 'r');

% Piirretään pystysuora viiva BPFO:n kohdalle
xline(BPFO, '--r', 'LineWidth', 1.5);
text(BPFO, 0.001*0.9, sprintf('BPFO', BPFO), ...
     'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
     'Color', 'r', 'FontSize', 9);

% Harmonisien viivojen värit
harm_color = [0.1, 0.6, 0.1]; 

% Taajuudet, joille piirretään harmoniset (voit säätää listaa)
base_freqs = [BPFO, BPFI, BSF, FTF, Fr];
labels = {'BPFO', 'BPFI', 'BSF', 'FTF', 'Fr'};



% Piirretään pystysuora viiva BPFI:n kohdalle
xline(BPFI, '--r', 'LineWidth', 1.5);
text(BPFI, 0.001*0.9, sprintf('BPFI', BPFI), ...
     'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
     'Color', 'r', 'FontSize', 9);

% Piirretään pystysuora viiva Fr:n kohdalle
xline(Fr, '--r', 'LineWidth', 1.5);
text(Fr, 0.001*0.9, sprintf('Fr', Fr), ...
     'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
     'Color', 'r', 'FontSize', 9);

% Piirretään pystysuora viiva BSF:n kohdalle
xline(BSF, '--r', 'LineWidth', 1.5);
text(BSF, 0.001*0.9, sprintf('BSF', BSF), ...
     'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
     'Color', 'r', 'FontSize', 9);

% Piirretään pystysuora viiva FTF:n kohdalle
xline(FTF, '--r', 'LineWidth', 1.5);
text(FTF, 0.001*0.9, sprintf('FTF', FTF), ...
     'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
     'Color', 'r', 'FontSize', 9);



% Harmoniset taajuudet

n_harmonics = 4;

for i = 1:length(base_freqs)
    for h = 2:n_harmonics  % Aloitetaan 2. harmonisesta
        f_harm = h * base_freqs(i);
        if f_harm <= max(f_plot)
            xline(f_harm, '--', 'Color', harm_color, 'LineWidth', 0.5);
            % Lisää teksti harmoniselle
            text(f_harm, -5, sprintf('%s %dx', labels{i}, h), ...
                 'Rotation', 90, 'VerticalAlignment', 'bottom', ...
                 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', harm_color);
        end
    end
end

title('FFT-spektri');
xlabel('Taajuus (Hz)');
ylabel('Amplitudi');
grid on;

%% Demoduloitu spektri (Envelope analysis)
% Band-pass suodatin (voit säätää rajataajuuksia)
bpFilt = designfilt('bandpassiir','FilterOrder',4, ...
         'HalfPowerFrequency1',500,'HalfPowerFrequency2',3000, ...
         'SampleRate',Fs);

% Suodatetaan signaali ja lasketaan envelope
signal_bp = filtfilt(bpFilt, signal);
env = abs(hilbert(signal_bp));

% FFT envelope-signaalista
N_env = length(env);
env_windowed = env .* hamming(N_env);
Y_env = fft(env_windowed);
f_env = Fs*(0:floor(N_env/2))/N_env;
P2_env = abs(Y_env/N_env);
P1_env = P2_env(1:floor(N_env/2)+1);
P1_env(2:end-1) = 2*P1_env(2:end-1);
P1_env_dB = 20*log10(P1_env);

% Rajataan taajuusalue
f_env_plot = f_env(f_env <= max_freq);
P1_env_plot_dB = P1_env_dB(f_env <= max_freq);

%% Piirretään molemmat spektrit samaan figureen


subplot(3,1,3);
plot(f_env_plot, P1_env_plot_dB, 'b', 'LineWidth', 1.2);
hold on;
ylim([-120, 0]);
title('Demoduloitu spektri (Envelope FFT)');
xlabel('Taajuus (Hz)');
ylabel('Amplitudi (dB)');
grid on;

% Taajuusviivat envelope-spektriin
for i = 1:length(base_freqs)
    xline(base_freqs(i), '--r', 'LineWidth', 1.2);
    text(base_freqs(i), -5, labels{i}, ...
         'Rotation', 90, 'VerticalAlignment', 'bottom', ...
         'HorizontalAlignment', 'right', 'Color', 'r', 'FontSize', 8);
    % Harmoniset komponentit (alkaen 2. harmonisesta)
for h = 2:n_harmonics
    f_harm = h * BPFO;
    if f_harm <= max(f_env_plot)
        xline(f_harm, '--', 'Color', harm_color, 'LineWidth', 0.5);
        text(f_harm, -10, sprintf('BPFO %dx', h), ...
            'Rotation', 90, 'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'center', 'FontSize', 8, ...
            'Color', harm_color);
    end
end
% Piirretään pystysuora viiva Fr:n ja BPFI:n modulaation kohdalle
xline(F1_mod_positive, '--r', 'LineWidth', 1.5);
text(F1_mod_positive, 0.001*0.9, sprintf('F1+', F1_mod_positive), ...
     'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
     'Color', 'r', 'FontSize', 9);
% Piirretään pystysuora viiva FTF:n ja BPFI:n modulaation kohdalle
xline(F2_mod_positive, '--r', 'LineWidth', 1.5);
text(F2_mod_positive, 0.001*0.9, sprintf('F2+', F2_mod_positive), ...
     'Rotation', 90, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
     'Color', 'r', 'FontSize', 9);


for h = 2:n_harmonics
    f_harm = h * BPFI;
    if f_harm <= max(f_env_plot)
        xline(f_harm, '--', 'Color', harm_color, 'LineWidth', 0.5);
        text(f_harm, -10, sprintf('BPFI %dx', h), ...
            'Rotation', 90, 'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'center', 'FontSize', 8, ...
            'Color', harm_color);
    end
end
for h = 2:n_harmonics
    f_harm = h * Fr;
    if f_harm <= max(f_env_plot)
        xline(f_harm, '--', 'Color', harm_color, 'LineWidth', 0.5);
        text(f_harm, -10, sprintf('Fr %dx', h), ...
            'Rotation', 90, 'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'center', 'FontSize', 8, ...
            'Color', harm_color);
    end
end
end

%% f/fs spektri 

% Normalisoitu taajuusakseli
f_norm = f_plot / Fr;

% Uusi kuvaaja
figure;
plot(f_norm, P1_dB, 'LineWidth', 1.5);
hold on;
plot(peak_freqs/Fr, 20*log10(peak_amps), 'ro', 'MarkerFaceColor', 'r');
grid on;
xlabel('f / f_r');
ylabel('Amplitudi (dB)');
title('FFT-spektri suhteutettuna pyörimisnopeuteen (f/fs)');

% Piirretään viivat tunnusomaisille taajuuksille
for i = 1:length(base_freqs)
    xline(base_freqs(i)/Fr, '--', labels{i}, ...
           'Color', harm_color, 'LineWidth', 1.2, ...
           'LabelOrientation', 'horizontal', ...
           'LabelVerticalAlignment', 'bottom');
end
