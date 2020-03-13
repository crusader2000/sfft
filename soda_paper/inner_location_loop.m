/*
 Step 1 of inner loop
 REQUIREMENTS:
 sigma and tau must lie in the range [0,n)
 sigma must be odd
*/
function [sigma,tau]=get_sigma_and_tau(n)
	tau=randi([0,n-1],1,1)
	sigma=2*randi([0,n/2-1],1,1)+1
end

/*
 Step 2 involves permuting and multiplying the signal with the window function

%INCOMPLETE PART
 Step 2 Part 1
 Creating a flat-window function of length n (the length of input signal)
*/
function g=get_window_function(n,epsilon_dash)
	g = chebwin(n);
	% BOX CAR FUNCTION
	box=rectangularPulse(-epsilon_dash*n,epsilon_dash*n,-n/2:n/2);
	%plot(-n/2:n/2,box)
	% Check the length of this
	if length(box)==length(g)
		g=cconv(box,g);
	else 
		display("length of box car function is not equal to length of the window function");
	end
end

% Step 2 part 2
% Multiplying the window function with the permuted signal
function g=permute(x,sigma,tau,n)
	% Returns a permuted signal
	% Assumption: The signal is periodic with period n
  i=0:n-1;
	g=x[(sigma*i+tau)%n];
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
	z=zeros(1,B);
	l=0:(w/B)-1;
	for i=1:B
		%check this line for dimesion errors
		%putting a mod might help
		tmp=i+B*(l);
		z[i] = sum(y[temp]);
	end
end

% Step 4: Making the hash function and offset function

function h_sigma=hash_function(n,sigma,B)
	h_sigma=round((sigma*B*(0:n-1))/n);
end

function o_sigma=offset_function(h_sigma,n,sigma,B):
	% h is the list returned by the hash function
	o_sigma=sigma*(0:n-1)-(n/B)*h_sigma;
end


% Step 5: The set which is the output in location loop
function I=finding_I(Z,n,h_sigma)
	max_element=max(Z);
	J=[];
	for i=0:length(Z)-1
		if Z[i]==max_element
			append(J,i);
		end
	end

	J=sort(J);
	I=[];

	for i=1:n
		[index]=binarySearch(J,length(J),h_sigma[i]);
		if index!=-1
			append(I,i-1);
		end
	end
end


%Step 6 : The estimates in estimation loop

function X_estim=find_estimates(h_sigma,o_sigma,tau,G,I,Z,n)
	X_estim=zeros(1,n);

	for j=1:length(I)
		i=I[j]+1;
		X_estim[i]=(Z[h_sigma[i]]*exp((1i*2*pi*tau*(i-1))/n))/G[o_sigma[i]];
	end

end


function [I,h_sigma,o_sigma,Z,tau,G]=inner_location_loop(x,n,k)

[sigma,tau]=get_sigma_and_tau(n,k);

%Make sure that B is able to divide n
B=sqrt(n*k);

epsilon = 1/B;

epsilon_dash = epsilon/2;

g=get_window_function(n,B,epsilon_dash)
G=fft(g)

x_permuted=permute(x,sigma,tau,n);
y=x_permuted.*g;

z=generate_z(y,w,B);
Z=fft(z);

h_sigma=hash_function(n,sigma,B);
o_sigma=offset_function(h_sigma,n,sigma,B);

I=finding_I(Z,n,h_sigma);

end