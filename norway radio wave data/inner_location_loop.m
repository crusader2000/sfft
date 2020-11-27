%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VVVIMP
% The values in I are in the range [0,n-1] both inclusive
%VVVIMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I,h_sigma,o_sigma,Z,tau] = inner_location_loop(x,n,k,g,d,B,delta)
    % clc;
    % close all;
    [sigma,tau]=get_sigma_and_tau(n);
    
    x_permuted=permute(x,sigma,tau,n);
    y=x_permuted.*g;

    w=round(0.7*B*log(n/delta));
    z=generate_z(y,w,B,n);
    Z=fft(z,B);

    h_sigma=hash_function(n,sigma,B);
    o_sigma=offset_function(h_sigma,n,sigma,B);

    I=finding_I(Z,n,d,k,h_sigma,B);
end

%  Step 1 of inner loop
%  REQUIREMENTS:
%  sigma and tau must lie in the range [0,n)
%  sigma must be 1 for the mapping criteria to be satisfied

function [sigma,tau]=get_sigma_and_tau(n)
    tau=randi([0,n-1],1,1);
    sigma= 2*randi([0,n/2-1],1,1)+1;
    % fprintf("  Randomly selected sigma and tau \n tau = %d sigma= %d\n",tau,sigma);
end

% Step 2
% Multiplying the window function with the permuted signal
function x_permute=permute(x,sigma,tau,n)
    % Returns a permuted signal
    % Assumption: The signal is periodic with period n
    ii=0:n-1;
    permute_order=int64(mod((sigma*ii+tau),n))+1;
    x_permute=x(permute_order);
end


%  Step 3 : Finding the DFT of the signal y with bucketisation
%  There are 2 steps involved here
%  Generating z
%  Finding the dft of z

function z=generate_z(y,w,B,n)
    % w is a parameter related to the window function
    % B is the bucket size
    
    z=zeros(1,B);
    ll=B*(0:round(w/B)-1);
    
    for k=0:B-1
        temp=k+ll;
        z(k+1)=sum(y(mod(temp,n)+1));
    end
end

% Step 4: Making the hash function and offset function
function h_sigma=hash_function(n,sigma,B)
    h_sigma=mod((0:sigma*B:sigma*B*(n-1)),B);
end

function o_sigma=offset_function(h_sigma,n,sigma,B)
    % h_sigma is the list returned by the hash function
    
    % o_sigma=mo1d((mod((sigma*(0:n-1)-((n/B)*h_sigma)),n)+n),n);

    h_sigma_dash=floor((0:sigma*B:sigma*B*(n-1))/n);
    o_sigma=sigma*(0:n-1)-((n/B)*h_sigma_dash);
end

% Step 5: The set which is the output in location loop
function I=finding_I(Z,n,d,k,h_sigma,B)

    if d*k > B
        disp("no of elements in J is more than the size of Z");
    else
        no_elements=d*k;    
    end

    [sort_Z,indices]=sort(abs(Z));
    
    % SORT IN DESCENDING ORDER OF MAGNITUDES 
    indices = flip(indices);

    J=indices(1:no_elements)-1;

    I=[];
    
    mask = zeros(1,n);

    for idx=1:length(J)
        temp = h_sigma == J(idx);
        mask = mask +temp;     
    end

    indices = 1:n;
    I = indices(mask>=1)-1;
end