% Step 1 run a location inner loop L times
% This function returns a set of Is (I1,I2,....,IL) and
% the sigma associated with it in the array set_sigma
function [set_I,set_sigma]=get_the_set_of_I(x,L)
	set_I=[]
	set_sigma=[]
	for i=1:L
		[I,sigma]=inner_location_loop(x,n)
		set_I[i]=I
		set_sigma[i]=sigma
	end
end

% Step 2 of outer loop
% We find the union of all Is in union_of_all_sets
% and find the required set s at the same time

function [s,union_of_all_sets]=get_s(I,r)	
	union_of_all_sets=[]
	for i=0:r-1
		I[i]=sort(I[i]])
		union_of_all_sets=union(union_of_all_sets,I[i])
	end
	s=zeros(max(union_of_all_sets),1)
	for i=1:len(union_of_all_sets)
		search_element=union_of_all_sets[i]
		for j=0:r-1
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
function X_r_I_dash=get_x_r_I_dash(I_dash,L)
	X_r_I_dash=[]
	for i=1:L
		X_estim=find_estimates(x,h_sigma,o_sigma,omega,tau,G,I_dash,Z)
	end
end