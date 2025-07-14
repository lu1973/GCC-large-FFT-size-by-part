[sig,fs] = audioread('data\226596.wav');
% create a shifted version of the input signal with 2 different delay values
refsig = [circshift(sig, -500);circshift(sig, -200)];
sig = [sig;sig];
audiowrite('data\226596_shifted.wav', refsig, fs);

buffer_length = 512;

[prm, states] = longGCC_init(fs, buffer_length, 42, 0.1, 300, 6000);

n_frame = round((length(sig)-buffer_length)/buffer_length);

delays = zeros(n_frame, 1);

for n=1:n_frame
    % the input signal is split into buffer of length 512 samples
    % but the GCC will be computed on a length of 21504 (42*512) samples
    sig_frm = sig(1+(n-1)*buffer_length:n*buffer_length);
    refsig_frm = refsig(1+(n-1)*buffer_length:n*buffer_length);
    [delay_out, states] = longGCC_process(sig_frm, refsig_frm, prm, states);
    delays(n) = states.delay;
endfor

figure()
t = linspace(0, (length(sig)-buffer_length)/fs, n_frame);
stairs(t, delays)
grid
xlabel('time (s)')
ylabel('delay (samples)')



