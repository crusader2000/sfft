function [I,h_sigma,o_sigma,Z,tau,G] = inner_location_loop(x,n,k)
    clc;
    close all;
    [sigma,tau]=get_sigma_and_tau(n);
    
    %Make sure that B is able to divide n
    B=sqrt(n*k);
    if mod(n,B) ~= 0
        B=B-mod(n,B);
    end
        
    epsilon = 1/B;
    epsilon_dash = epsilon/2;
    
    g=get_window_function(n,epsilon_dash);
    G=fft(g);
    
    x_permuted=permute(x,sigma,tau,n);
    y=x_permuted.*g;
%     plot(1:n,y,1:n,x,'o');
%     legend('y','x');
    % w = O(B log( n/delta)).
    % delta= 1/n^c. choose c so that n is negligibly small.
    w=B*(n+20);
    z=generate_z(y,w,B,n);
    Z=fft(z);

    h_sigma=hash_function(n,sigma,B);
    o_sigma=offset_function(h_sigma,n,sigma,B);
%     disp("h_sigma");
%     disp(h_sigma);
%     disp("o_sigma");
%     disp(o_sigma);
    
    %since d is in the order of 1/epsilon
    d=int64(1/epsilon);
    
%     disp(abs(Z));
%     disp(length(Z));
%     disp(d*k*n/B)    
    % Problem conjugate numbers hav the same magnitude. This will
    % interfere with elements in J and finding them 
%    I=zeros(1,d*k);
    I=finding_I(Z,n,d,k,h_sigma);
end

%  Step 1 of inner loop
%  REQUIREMENTS:
%  sigma and tau must lie in the range [0,n)
%  sigma must be odd

function [sigma,tau]=get_sigma_and_tau(n)
    tau=randi([0,n-1],1,1);
    sigma=2*randi([0,n/2-1],1,1)+1;
    fprintf("  Randomly selected sigma and tau \n tau = %d sigma= %d\n",tau,sigma);
end


%  Step 2 involves permuting and multiplying the signal with the window function
%
%   Step 2 Part 1
%  Creating a flat-window function of length n (the length of input signal)

function g=get_window_function(n,epsilon_dash)
    chewin = chebwin(n);
    % BOX CAR FUNCTION
    box=rectangularPulse(-epsilon_dash*n,epsilon_dash*n,((-n/2):(n/2)-1));
    g=[];
    % Check the length of this
    if length(box)==length(chewin)
        g=cconv(box,chewin,n);
%         disp("Plotting the window function")
%         subplot(3,1,1);
%         plot(-n/2:n/2-1,box);
%         title('box car function')
%         subplot(3,1,2);
%         plot(-n/2:n/2-1,chewin);
%         title('dolph-chebyshev function')
%         subplot(3,1,3);
%         plot(-n/2:n/2-1,g);
%         title('standard window function')
%         figure;
    else
        disp("length of box car function is not equal to length of the chebwin window function");
    end
end

% Step 2 part 2
% Multiplying the window function with the permuted signal
function x_permute=permute(x,sigma,tau,n)
    % Returns a permuted signal
    % Assumption: The signal is periodic with period n
    i=0:n-1;
    temp=int64(mod((sigma*i+tau),n));
    %        disp(temp);
    %     disp(length(temp));
    %     disp(length(x))
    x_permute=zeros(1,n);
    for k = 1:n
        x_permute(k) = x(temp(k)+1);
        %display(k+" "+length(x_permute))
        %display(mod(temp(k),n))
    end
%     disp("Plotting the permuted signal");
%     stem(1:n,x);
%     hold on;
%     stem(1:n,x_permute);
%     legend('x','x_permute');
%     hold off;
%     figure;
end


%  Step 3 : Finding the DFT of the signal y with bucketisation
%  There are 2 steps involved here
%  Generating z
%  Finding the dft of z

function z=generate_z(y,w,B,n)
    % w is a parameter related to the window function
    % B is the bucket size
    z=zeros(1,B);
    l=B*(0:(w/B)-1);
%     disp(l);
    for i=0:B-1
        %check this line for dimesion errors
        %putting a mod might help
        temp=i+l;
        z(i+1)=sum(y(mod(temp,n)+1));
    end
%     stem(1:B,z);
%     title('z in time domain');
%     figure;
end

% Step 4: Making the hash function and offset function
function h_sigma=hash_function(n,sigma,B)
    h_sigma=round((sigma*B*(0:n-1))/n);
%     plot(0:n-1,h_sigma);
%     figure;
end

function o_sigma=offset_function(h_sigma,n,sigma,B)
    % h is the list returned by the hash function
    o_sigma=sigma*(0:n-1)-(n/B)*h_sigma;
%     plot(0:n-1,o_sigma);
%     figure;
end

% Step 5: The set which is the output in location loop
function I=finding_I(Z,n,d,k,h_sigma)
    J=[];
    sorted_Z=sort(abs(Z));
    
    
%     no_elements=int64(d*k);
     no_elements=4;    
%     uncomment the above line once d and k are properly assigned


    % Taking only the sqrt(dk) terms into consideration
    [sort_Z,indices]=sort(Z);
    J=indices(1:no_elements);
%     Z(J)
    J=sort(J);
    I=[];
    
    for i=1:n
        [index]=binarySearch(J,length(J),h_sigma(i));
        if index ~= -1
            I=[I,i-1];
        end
    end
end
