function X = fft_recur(x)
    n = length(x);
    X = zeros(n,1);
    if n == 2
        X(1)= x(1)+x(2);
        X(2)= x(1)-x(2);
    else
    X_even = fft_recur(x(1:2:n));
    X_odd = fft_recur(x(2:2:n));
    W =  exp(-2i*pi*(0:(n/2-1))/n)';
    X = [X_even+W.*X_odd;X_even-W.*X_odd];
    end
end