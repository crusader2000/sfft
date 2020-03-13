% Step 1 run a location inner loop L times
% This function returns a set of Is (I1,I2,....,IL) and
% the sigma associated with it in the array set_sigma
function [set_I,set_h_sigma,set_o_sigma,set_Z,set_tau,set_G]=get_the_sets(x,L,n)
	set_I=[]
	set_h_sigma=[]
	set_o_sigma=[]
	set_tau=[]
	set_Z=[]
	set_G=[]
	for i=1:L
		[I,h_sigma,o_sigma,Z,tau,G]=inner_location_loop(x,n)
		append(set_I,I)
		append(set_h_sigma,h_sigma)
		append(set_o_sigma,o_sigma)
		append(set_Z,Z)
		append(set_tau,tau)
		append(set_G,G)
	end
end

% Step 2 of outer loop
% We find the union of all Is in union_of_all_sets
% and find the required set s at the same time

function [s,union_of_all_sets]=get_s(I,L)
	union_of_all_sets=[]
	for i=1:L
		I[i]=sort(I[i])
		union_of_all_sets=union(union_of_all_sets,I[i])
	end
	s=zeros(max(union_of_all_sets),1)
	for i=1:len(union_of_all_sets)
		search_element=union_of_all_sets[i]
		for j=1:L
			[index]=binarySearch(I[j],length(I[j]),search_element)
			if index!=-1
				s[search_element]=s[search_element]+1
			end
		end
	end
end

% Step 3 of outer loop
% We find I_dash as per the formula given
function I_dash=get_I_dash(union_of_all_sets,s,L)
	I_dash=[]
	for i=1:length(union_of_all_sets)
		if s[union_of_all_sets[i]]>= L/2
			append(I_dash,union_of_all_sets[i])
		end
	end
end

% Step 4 of outer loop
% We run the estimation inner loop L times to get X_r_I_dash
function set_of_X_r_I_dash=get_x_r_I_dash(I_dash,set_h_sigma,set_o_sigma,set_Z,set_tau,set_G)
	set_of_X_r_I_dash=[]
	for i=1:L
		X_estim=find_estimates(set_h_sigma[i],set_o_sigma[i],set_tau[i],set_G[i],I_dash,set_Z[i])
		append(set_of_X_r_I_dash,X_estim)
	end
end


% Step 5 of outer loop
% We need to find the median of real and imaginary parts

function X=get_the_fourier_transform(set_of_X_r_I_dash,I_dash)
	X=zeros(length(set_of_X_r_I_dash[1]),1)
	for i=1:length(I_dash)
		X_i_real=[]
		X_i_imag=[]
		for j=1:size(set_of_X_r_I_dash,1)
			append(X_i_real,real(set_of_X_r_I_dash[j][I_dash[i]])
			append(X_i_imag,imag(set_of_X_r_I_dash[j][I_dash[i]])
		end
		X[i]=median(X_i_real)+j*median(X_i_imag)
	end
end

function X=sfft(x,L,n)
		[set_I,set_h_sigma,set_o_sigma,set_Z,set_tau,set_G]=get_the_sets(x,L,n)
		[s,union_of_all_sets]=get_s(I,L)
		I_dash=get_I_dash(union_of_all_sets,s,L)
		set_of_X_r_I_dash=get_x_r_I_dash(I_dash,set_h_sigma,set_o_sigma,set_Z,set_tau,set_G)
		X=get_the_fourier_transform(set_of_X_r_I_dash,I_dash)
end
