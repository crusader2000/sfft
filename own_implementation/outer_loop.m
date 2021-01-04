function X=outer_loop(x,n,k)
    %since L is in the order of log(n)
    L=int64(log2(n));
    
    [set_I,set_h_sigma,set_o_sigma,set_Z,set_tau,set_G]=get_the_sets(x,L,n,k);
    close all;
    [s,union_of_all_sets]=get_s(set_I,L);
    I_dash=get_I_dash(union_of_all_sets,s,L)
    disp(set_I)
    set_of_X_r_I_dash=get_x_r_I_dash(I_dash,set_h_sigma,set_o_sigma,set_Z,set_tau,set_G,L,n);
%       disp(set_of_X_r_I_dash)
     X=get_the_fourier_transform(set_of_X_r_I_dash,I_dash,n);
     plot(X);
    
end

% Step 1 run a location inner loop L times
% This function returns a set of Is (I1,I2,....,IL) and
% the sigma associated with it in the array set_sigma
function [set_I,set_h_sigma,set_o_sigma,set_Z,set_tau,set_G]=get_the_sets(x,L,n,k)
    set_I=[];
    set_h_sigma=[];
    set_o_sigma=[];
    set_tau=[];
    set_Z=[];
    set_G=[];
    for i=1:L
        [I,h_sigma,o_sigma,Z,tau,G]=inner_location_loop(x,n,k);
        set_I{i}=I;
        set_h_sigma=[set_h_sigma;h_sigma];
        set_o_sigma=[set_o_sigma;o_sigma];
        set_Z=[set_Z;Z];
        set_tau=[set_tau;tau];
        set_G=[set_G;G];
    end
    disp("Exiting the get_the_sets function");
%     disp(set_h_sigma(1,:));
end

% Step 2 of outer loop
% We find the union of all Is in union_of_all_sets
% and find the required set s at the same time

function [s,union_of_all_sets]=get_s(I,L)
    union_of_all_sets=[];
    for i=1:L
%         disp(I{i})
        union_of_all_sets=union(union_of_all_sets,I{i});
    end
    disp("Found the union of all sets");
    
    % Need to include zero too
    s=zeros(max(union_of_all_sets)+1,1);
    for r=1:L
        I_r=I{r};
        elements=length(I_r);
        for j=1:elements
                s(I_r(j)+1)=s(I_r(j)+1)+1;
        end
    end
    disp("Found s");
end

% Step 3 of outer loop
% We find I_dash as per the formula given
function I_dash=get_I_dash(union_of_all_sets,s,L)
    I_dash=[];
    for i=1:length(union_of_all_sets)
        if s(union_of_all_sets(i)+1)>= L/2
            I_dash=[I_dash;union_of_all_sets(i)];
        end
    end
   disp("Found I_dash");
end

% Step 4 of outer loop
% We run the estimation inner loop L times to get X_r_I_dash
function set_of_X_r_I_dash=get_x_r_I_dash(I_dash,set_h_sigma,set_o_sigma,set_Z,set_tau,set_G,L,n)
    set_of_X_r_I_dash=[];
    for i=1:L
        X_estim=find_estimates(set_h_sigma(i,:),set_o_sigma(i,:),set_tau(i),set_G(i,:),I_dash,set_Z(i,:),n);
%         disp(X_estim);
        set_of_X_r_I_dash{i}=X_estim;
    end
    disp("found set_of_X_r_I_dash");
   
end

% Step 5 of outer loop
% We need to find the median of real and imaginary parts

function X=get_the_fourier_transform(set_of_X_r_I_dash,I_dash,n)
    X=zeros(n,1);
    for i=1:length(I_dash)
        X_i_real=set_of_X_r_I_dash(:,i);
        X_i_imag=set_of_X_r_I_dash(:,i);
        X(i)=median(X_i_real)+1i*median(X_i_imag);
    end
end


function X_estim=find_estimates(h_sigma,o_sigma,tau,G,I,Z,n)
    X_estim=zeros(1,n);
%     disp("in find_sdtimates")
%     length(I)
    for j=1:length(I)
        k=I(j);
        k=mod(k,n)+1;
        X_estim(k)=(Z(h_sigma(k))*exp((1i*2*pi*tau*(k-1))/n))/G(o_sigma(k));
%         disp(X_estim(k))
    end
end
