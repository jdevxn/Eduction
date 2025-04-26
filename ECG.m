% ECG Signal Generation and Analysis
clear;
close all;
clc;

% Time parameters
fs = 1000; % Sampling frequency (Hz)
T = 10; % Duration (s)

% Generate dynamic ECG signal
[ecg, t] = generate_dynamic_ecg(fs, T);

% Add baseline wander and noise to the ECG signal
ecg_noisy = add_noise(ecg, t);

% Filter noisy ECG signal using a bandpass Butterworth filter
ecg_filtered = filter_ecg(ecg_noisy, fs);

% Detect R-peaks and calculate heart rate
[hr_bpm, r_peaks, rr_intervals] = calculate_heart_rate(ecg_filtered, fs);

% AFib Detection based on RR interval irregularity
if isempty(rr_intervals) || any(isnan(rr_intervals))
    disp('Error: Invalid RR intervals.');
else
    rr_diff = abs(diff(rr_intervals)); % Difference between successive RR intervals
    rr_std = std(rr_diff); % Standard deviation of RR interval differences
    fprintf('Standard deviation of RR intervals: %.4f seconds\n', rr_std);

    % Define AFib threshold
    afib_threshold = 0.12; % 120 ms

    % Check for AFib
    if rr_std > afib_threshold
        disp('Warning: Possible AFib detected!');
    else
        disp('Normal rhythm detected.');
    end

    % Display average heart rate
    fprintf('Average heart rate: %.2f bpm\n', mean(hr_bpm));
end

% Detect R-peaks and calculate heart rate
function [hr_bpm, r_peaks, rr_intervals] = calculate_heart_rate(ecg_filtered, fs)
    threshold = 0.7; % R-peak detection threshold
    [~, r_peaks] = findpeaks(ecg_filtered, 'MinPeakHeight', threshold);
    
    if isempty(r_peaks)
        % Handle case when no R-peaks are detected
        warning('No R-peaks detected. Heart rate cannot be calculated.');
        hr_bpm = NaN; % Return NaN if no R-peaks detected
        rr_intervals = NaN;
        return;
    end
    
    hr_samples = diff(r_peaks); % Difference between R-peak locations in samples
    rr_intervals = hr_samples / fs; % R-R intervals (seconds)
    hr_bpm = 60 ./ rr_intervals; % Heart rate (bpm)
end

% Plot the ECG signals
plot_ecg_signals(t, ecg_noisy, ecg_filtered, hr_bpm, r_peaks, fs);

% Generate dynamic ECG signal
function [ecg, t] = generate_dynamic_ecg(fs, T)
    hr_min = 60; % Minimum heart rate (bpm)
    hr_max = 100; % Maximum heart rate (bpm)
    num_cycles = T*hr_min/60; % Calculate the approximate number of cycles
    t_hr = linspace(0, T, num_cycles); % Time vector for heart rate changes
    hr = interp1([0, T], [hr_min, hr_max], t_hr, 'pchip'); % Interpolate between min and max heart rate

    % Preallocate ECG array
    ecg = zeros(1, fs * T);

    idx = 1;
    for i = 1:length(hr)
        f_ecg = hr(i) / 60; % ECG frequency (Hz)
        ecg_cycle = round(fs/f_ecg); % ECG cycle length in samples
        ecg_temp = ecg_cycle_gen(linspace(0, 1, ecg_cycle)); % Generate ECG cycle
        ecg(idx:idx + length(ecg_temp) - 1) = ecg_temp; % Append ECG cycle to signal
        idx = idx + length(ecg_temp);
    end

    % Trim ECG array to remove excess zeros
    ecg = ecg(1:find(ecg, 1, 'last'));

    t = 0:1/fs:(length(ecg)-1)/fs; % Adjust the time vector to match the ECG signal length
end

% Add baseline wander and noise to the ECG signal
function ecg_noisy = add_noise(ecg, t)
    bw_freq = 0.1; % Baseline wander frequency (Hz)
    bw_amp = 0.05; % Baseline wander amplitude
    baseline_wander = bw_amp * sin(2*pi*bw_freq*t);
    noise_amp = 0.05; % Noise amplitude
    noise = noise_amp * randn(1, length(t));
    ecg_noisy = ecg + baseline_wander + noise;
end

% Filter noisy ECG signal using a bandpass Butterworth filter
function ecg_filtered = filter_ecg(ecg_noisy, fs)
    order = 2; % Filter order
    low_cutoff = 0.5; % Low cutoff frequency (Hz)
    high_cutoff = 50; % High cutoff frequency (Hz)
    Wn = [low_cutoff, high_cutoff] / (fs/2); % Normalized cutoff frequencies
    [b, a] = butter(order, Wn, 'bandpass'); % Design Butterworth filter
    ecg_filtered = filter(b, a, ecg_noisy);
end

% Plot the ECG signals
function plot_ecg_signals(t, ecg_noisy, ecg_filtered, hr_bpm, r_peaks, fs)
    figure;
    subplot(3, 1, 1);
    plot(t, ecg_noisy);
    title('Noisy ECG Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');

    subplot(3, 1, 2);
    plot(t, ecg_filtered);
    hold on;
    plot(r_peaks/fs, ecg_filtered(r_peaks), 'ro');
    title('Filtered ECG Signal and R-Peaks');
    xlabel('Time (s)');
    ylabel('Amplitude');

    subplot(3, 1, 3);
    plot(t(r_peaks(1:end-1)), hr_bpm);
    title('Heart Rate');
    xlabel('Time (s)');
    ylabel('BPM');
end

% Generate a single ECG cycle
function y = ecg_cycle_gen(t)
    % P wave
    p_wave_amp = 0.2;
    p_wave_dur = 0.08;
    p_wave_t = t - 0.2;
    p_wave = p_wave_amp * exp(-(p_wave_t/p_wave_dur).^2);
    
    % QRS complex
    qrs_amp = 1;
    qrs_dur = 0.1;
    qrs_t = t - 0.4;
    qrs = qrs_amp * exp(-(qrs_t/qrs_dur).^2);
    
    % T wave
    t_wave_amp = 0.3;
    t_wave_dur = 0.12;
    t_wave_t = t - 0.6;
    t_wave = t_wave_amp * exp(-(t_wave_t/t_wave_dur).^2);
    
    % Combine P, QRS, and T waves
    y = p_wave + qrs + t_wave;
end
