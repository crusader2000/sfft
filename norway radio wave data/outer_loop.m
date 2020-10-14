function X=outer_loop(x,n,k,g,d,B,delta,G)
    %since L is in the order of log(n)
    L = int64(log2(n));
   
    [set_I,set_h_sigma,set_o_sigma,set_Z,set_tau] = get_the_sets(x,L,n,k,g,d,B,delta);
  
    [s,union_of_all_sets] = get_s(set_I,L);
  
    I_dash = get_I_dash(union_of_all_sets,s,L);
  
    set_of_X_r_I_dash = get_x_r_I_dash(I_dash,set_h_sigma,set_o_sigma,set_Z,set_tau,G,L,n);
  
    X = get_the_fourier_transform(set_of_X_r_I_dash,I_dash,n);
    % plot(1:n,X);
    % figure;
    % dBFS_X = 20*log10(abs(X)) - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;
    % plot(1:n,dBFS_X);
    
end

% Step 1 run a location inner loop L times
% This function returns a set of Is (I1,I2,....,IL) and
% the sigma associated with it in the array set_sigma
function [set_I,set_h_sigma,set_o_sigma,set_Z,set_tau]=get_the_sets(x,L,n,k,g,d,B,delta)
    set_I=[];
    set_h_sigma=[];
    set_o_sigma=[];
    set_tau=[];
    set_Z=[];

    for i=1:L
        [I,h_sigma,o_sigma,Z,tau]=inner_location_loop(x,n,k,g,d,B,delta);
        set_I=[set_I ;I];
        set_h_sigma=[set_h_sigma;h_sigma];
        set_o_sigma=[set_o_sigma;o_sigma];
        set_Z=[set_Z;Z];
        set_tau=[set_tau;tau];
    end
end

% Step 2 of outer loop
% We find the union of all Is in union_of_all_sets
% and find the required set s at the same time

function [s,union_of_all_sets]=get_s(I,L)
    union_of_all_sets=[];
    for i=1:L
        union_of_all_sets=union(union_of_all_sets,I(i,:));
    end

    s=zeros(max(union_of_all_sets)+1,1);
    for r=1:L
        I_r=I(r,:);
        num_elements=length(I_r);
        s(I_r(1:num_elements)+1)=s(I_r(1:num_elements)+1)+1;
    end
end

% Step 3 of outer loop
% We find I_dash as per the formula given
function I_dash=get_I_dash(union_of_all_sets,s,L)
    I_dash=[];
 
    for idx=1:length(union_of_all_sets)
        if s(union_of_all_sets(idx)+1)>= L/2
            I_dash=[I_dash;union_of_all_sets(idx)];
        end
    end
end

% Step 4 of outer loop
% We run the estimation inner loop L times to get X_r_I_dash
function set_of_X_r_I_dash=get_x_r_I_dash(I_dash,set_h_sigma,set_o_sigma,set_Z,set_tau,G,L,n)
    set_of_X_r_I_dash=[];
    for i=1:L
        X_estim=find_estimates(set_h_sigma(i,:),set_o_sigma(i,:),set_tau(i),G,I_dash,set_Z(i,:),n);
        set_of_X_r_I_dash=[set_of_X_r_I_dash;X_estim];
    end
end

% Step 5 of outer loop
% We need to find the median of real and imaginary parts

function X=get_the_fourier_transform(set_of_X_r_I_dash,I_dash,n)
    X=zeros(n,1);
    for k=1:length(I_dash)
        X_i_real=real(set_of_X_r_I_dash(:,int64(I_dash(k)+1)));
        X_i_imag=imag(set_of_X_r_I_dash(:,int64(I_dash(k)+1)));
        X(k)=median(X_i_real)+1i*median(X_i_imag);
    end

end


function X_estim=find_estimates(h_sigma,o_sigma,tau,G,I,Z,n)
    X_estim=zeros(1,n);

    for idx=1:length(I)
        k=I(idx)+1;
        X_estim(k)=(Z(h_sigma(k)+1)*exp((1i*2*pi*tau*(k-1))/n))/G(o_sigma(k)+1);
    end

end
