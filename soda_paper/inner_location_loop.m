/*
 Step 1 of inner loop
 REQUIREMENTS:
 sigma and tau must lie in the range [0,n)
 sigma must be odd
*/
function [sigma,tau]=get_sigma_and_tau(n)
	tau=randi([0,n-2],1,1)
	sigma=randi([0,n/2-1],1,1)
end

/*
 Step 2 involves permuting and multiplying the signal with the window function

 Step 2 Part 1
 Creating a flat-window function of length n (the length of input signal)
*/

% Step 2 part 2
% Multiplying the window function with the permuted signal
function g=permute(x,sigma,tau)
	% Returns a permuted signal
	% Assumption: The signal is periodic with period n
	n=length(x)
	g=x[(sigma*i+tau)%n]
end




function  y=multiplication_with_window_function(x,g)
	% x is the signal 
	% g is the window function
	y=x.*g
end

/*
# Step 3 : Finding the DFT of the signal y with bucketisation
# There are 2 steps involved here 
# Generating z
# Finding the dft of z 
*/

function z=generate_z(y,w,B)
	% w is a parameter related to the window function
	% B is the bucket size
	z=zeros(1,n)
	l=0:(w/B)-1
	for i=0:B 
		%check this line for dimesion errors
		tmp=i+B*(l)
		z[i] = sum(y[temp])
end


function Z=compute_dft(z)
	Z=fft(Z)
end


% Step 4: Making the hash function and offset function

function h_sigma=hash_function(n,sigma,B)
	h_sigma=round((sigma*B*(0:n-1))/n)
end

function o_sigma=offset_function(h_sigma,n,sigma,B):
	% h is the list returned by the hash function
	o_sigma=sigma*(0:n-1)-(n/B)*h_sigma
end

function I=finding_J(Z,n)
	max_element=max(Z)
	J=[]
	for i=0:length(Z)-1
		if Z[i]==max_element
			append(J,i)
		end
	end

%Incomplete function



end