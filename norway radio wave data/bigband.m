function X = bigband(x,coprimes,n)
    close all;
    X = zeros(1,n);

    num_coprimes = length(coprimes);
    
    b = {};
    b_dash = {};

    for ii = 1:num_coprimes
        indices = 0:coprimes(ii):n-1;
        b{ii} = fft(x(indices+1));
     
        b_dash{ii} = fft(x(indices+2));
        % plot(1:length(b{ii}),abs(b{ii}));
        % figure;
        % plot(1:length(b_dash{ii}),abs(b_dash{ii}));
        % figure;
    end

    num_iters = log2(20);
    
    for ii = 1:num_iters
        for jj = 1:num_coprimes
            for kk = 1:length(b{jj})
                % disp(kk);
                if abs(b{jj}(kk)) == abs(b_dash{jj}(kk))
                    f = ((angle(b_dash{jj}(kk)) - angle(b{jj}(kk)))*n)/(2*pi);
                    f = mod(f,n);
                    disp("f");
                    disp(f);
                   
                    X(int64(f)+1) = b{jj}(kk);

                    for m = 1:num_coprimes
                        f_dash = mod(int64(f),n/coprimes(m));
                        disp("f_dash");
                        disp(f_dash);
                        b{m}(f_dash+1) = b{m}(f_dash+1) - X(int64(f)+1);
                        b_dash{m}(f_dash+1) = b_dash{m}(f_dash+1) - X(int64(f)+1)*exp((1i*2*pi*f)/n);
                        % plot(1:length(b{m}),abs(b{m}));
                        % title('b{m}')
                        % figure;
                        % plot(1:length(b_dash{m}),abs(b_dash{m}));
                        % title('b_dash{m}')
                        % figure;
                    end
                end
            end
        end
    end
    % temp = abs(fft(x));
    % disp(temp(temp~=0));
    plot(1:n,abs(fft(x)))
    hold on;
    plot(1:n,abs(X));
    legend('Inbuilt','BigBand');
end