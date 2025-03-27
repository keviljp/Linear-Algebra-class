% define the signals and their periods
signals = {
    @(n) 3 + sin((4*pi/5)*n + pi/10) + cos(2*pi*n) + (-1).^n, 10; 
    @(n) sum_kronecker_delta(n, 5) - 2*sum_kronecker_delta(n-2, 5), 5;
    @(n) 1 - sin((pi/2)*n), 4; 
    @(n) 1 - sin((pi/2)*n), 16 
};

for idx = 1:size(signals, 1)
    %do fft for each and normalize with 1/N for DTFS
    x_func = signals{idx, 1};
    N = signals{idx, 2};
    
    n = 0:N-1;
    x = x_func(n);
    
    Cm = (1/N) * fft(x);
    
    % plot results
    figure;
    
    % plot x[n]
    subplot(3, 1, 1);
    stem(n, x, 'filled', 'LineWidth', 1.5);
    xlabel('n');
    ylabel('x[n]');
    title(['Signal x[n] (N = ', num2str(N), ')']);
    grid on;
    
    % plot |Cm|
    subplot(3, 1, 2);
    stem(n, abs(Cm), 'filled', 'LineWidth', 1.5);
    xlabel('m');
    ylabel('|C_m|');
    title(['Magnitude of C_m (N = ', num2str(N), ')']);
    grid on;
    
    % plot angle of Cm
    subplot(3, 1, 3);
    stem(n, angle(Cm), 'filled', 'LineWidth', 1.5);
    xlabel('m');
    ylabel('Angle(C_m)');
    title(['Phase of C_m (N = ', num2str(N), ')']);
    grid on;
end

% function for kronecker delta
function x = sum_kronecker_delta(n, K)
    x = zeros(size(n));
    for k = -100:100 % big range to approximate all of Z
        x = x + (n == k*K) - 2*(n == 2 + k*K);
    end
end