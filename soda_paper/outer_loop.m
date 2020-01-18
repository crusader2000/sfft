% Step 2 of outer loop
function [s,union_of_all_sets]=get_s(I,r)
	
	union_of_all_sets=[]
	for i=0:r-1
		I[i]=sort(I[i]])
		union_of_all_sets=union(union_of_all_sets,I[i])
	end
	s=zeros(max(union_of_all_sets)+1,1)
	for i=0:len(union_of_all_sets)-1
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
function I_dash=get_I_dash(I,union_of_all_sets,s,L)
	I_dash=[]
	for i=0:len(union_of_all_sets)-1
		if s[union_of_all_sets[i]]>= L/2
			append(I_dash,union_of_all_sets[i])
		end
	end
end
