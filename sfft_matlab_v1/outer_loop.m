function oloop_output = outer_loop(origx,n, filter_timedo,filter_sizet,filter_freqdo,filter_est_timedo,filter_est_sizet,filter_est_freqdo, B2,num,B,W_Comb,Comb_loops,loop_threshold,location_loops, loops)
  % disp("STARTING OUTER LOOP");
  permute_vals = zeros(1,loops);
  permuteb_vals = zeros(1,loops);
  disp(B);
  disp(B2);
  x_samp = {};
  for ii = 0:loops
    if ii<location_loops
    x_samp{ii+1} = zeros(1,B);
    else
    x_samp{ii+1} = zeros(1,B2);
    end
  end
 
  hits_found = 0;

  score = zeros(1,n);
  hits = zeros(1,n);
  
  %  BEGIN INNER LOOPS
  for ii=0:loops-1
    % disp("RUNNING INNER LOOP NUMBER");
    % disp(ii+1);

    a = 0;
    % b = mod(randi(n),n);
    b = 0;

    while gcd(a, n) ~= 1

      a = int32(mod(randi(n),n));
    end

    ai = mod_inverse(a, n);
    % fprintf("ai=%d , a=%d, n=%d, gcd=%d, mod=%d\n",ai,a,n,gcd(a,n),mod(ai*a,n));

    permute_vals(ii+1) = ai;
    permuteb_vals(ii+1) = b;

    % perform_location = ii < location_loops;

    % if perform_location
    %   cur_B = B;
    %   cur_filter_timedo = filter_timedo;
    %   cur_filter_sizet = filter_sizet;
    %   cur_filter_freqdo = filter_freqdo;
    % else
      cur_B = B2;
      cur_filter_timedo = filter_est_timedo;
      cur_filter_sizet = filter_est_sizet;
      cur_filter_freqdo = filter_est_freqdo;  
    % end

    J = zeros(1,num);

    [x_samp{ii+1},J] = inner_loop_locate(origx, n, cur_filter_timedo,cur_filter_sizet,cur_filter_freqdo, num, cur_B,a, ai, b);

    % if perform_location 
      [score, hits, hits_found] = inner_loop_filter_regular(J, n, num, cur_B,a, ai, b, loop_threshold,score, hits, hits_found);
      
    % end
  end
  
  % disp("unique indices");
  % disp(sum~=0));
  % disp("");

  
  % disp("loop threshold")
  % disp(loop_threshold);
  % disp("hits");
  % disp(hits);
  % disp("hits_found");
  % disp(hits_found);

  
  %END INNER LOOPS
  oloop_output = estimate_values(hits,hits_found,x_samp, loops,n, permute_vals, permuteb_vals, B, B2,filter_timedo,filter_sizet,filter_freqdo, filter_est_timedo,filter_est_sizet,filter_est_freqdo,location_loops);
end

function [x_samp,J] = inner_loop_locate(origx,n, filter_timedo,filter_sizet,filter_freqdo,num,B,a,ai,b)
  if mod(n,B) ~= 0
    disp("Warning: n is not divisible by B, which algorithm expects.\n");
  end
  
  x_sampt = zeros(1,B);
  % disp(filter_sizet);
  index=b;
  temp = zeros(1,filter_sizet);
  % fprintf("Inner loop locate\n");
  for ii = 0:filter_sizet-1
    % fprintf("index %d origx value %d\n",index,origx(index+1));
    x_sampt(mod(ii,B)+1) = x_sampt(mod(ii,B)+1) + origx(index+1)*filter_timedo(ii+1);
    index = mod((index+ai),n);
  end
  % disp(temp);
  % fprintf("\n");

  % disp("INNER LOOP LOCATE FFT_RECUR");
  x_samp = fft_recur(x_sampt);
  % disp(x_samp);
  % disp(size(x_samp));
  % disp(B);
  samples = abs(x_samp(1:B));
  
  % disp("INNER LOOP LOCATE FIND_LARGEST_INDICES");

  if num < B
    J= find_largest_indices(num,samples);
  else
    J= find_largest_indices(B,samples);
  end
end

function [score, hits, hits_found] = inner_loop_filter_regular(J,n,num,B,a,ai,b,loop_threshold,score, hits, hits_found)
  for ii=0:length(J)-1
    low = mod((int32(ceil((J(ii+1) - 0.5) * n / B)) + n),n);
    high = mod((int32(ceil((J(ii+1) + 0.5) * n / B)) + n),n);
    loc = timesmod(low, a, n);

    % fprintf("J(ii+1) = %d, low = %d, high = %d, loc = %d\n",J(ii+1),low,high,loc);
    
    jj= low;
    while jj ~= high
      score(loc+1) = score(loc+1)+1;

        if score(loc+1)==loop_threshold 
          hits(hits_found+1)=loc;
          hits_found = hits_found + 1;
        end
        loc = mod((loc + a),n);

      jj = mod((jj+1),n);
    end
  end
end

function oloop_output = estimate_values(hits,hits_found,x_samp, loops,n, permute_vals, permuteb_vals, B, B2,filter_timedo,filter_sizet,filter_freqdo, filter_est_timedo,filter_est_sizet,filter_est_freqdo,location_loops)
  
  oloop_output = zeros(1,n);
  values = zeros(2,loops);

  for ii = 0:hits_found-1
    position = 0;

    for jj = 0:loops-1
      if jj <location_loops
        cur_B = B;
        cur_filter_timedo = filter_timedo;
        cur_filter_sizet = filter_sizet;
        cur_filter_freqdo = filter_freqdo;
      else
        cur_B = B2;
        cur_filter_timedo = filter_est_timedo;
        cur_filter_sizet = filter_est_sizet;
        cur_filter_freqdo = filter_est_freqdo;  
      end
      permuted_index= timesmod(permute_vals(jj+1), hits(ii+1),  n);
      permuted_index = int32(permuted_index);
      cur_B = int32(cur_B);
      % n = int32(n);
     
      hashed_to = ((permuted_index*cur_B)/int32(n));
      dista = mod(permuted_index , (int32(n) / cur_B));
     
      if dista > ((n/cur_B)/2) 
        hashed_to = hashed_to + 1;
        dista = (dista -(n/cur_B));
      end
      dista = mod((n - dista),n);

      hashed_to = mod(hashed_to ,cur_B);
      filter_value = cur_filter_freqdo(dista+1);

      % exp_factor = exp((1i*2*pi*timesmod(permuteb_vals(jj+1), hits(ii+1), n)/n));
      exp_factor = 1;
      values(1,position+1) = real((x_samp{jj+1}(hashed_to+1)*exp_factor) / filter_value);
      values(2,position+1) = imag((x_samp{jj+1}(hashed_to+1)*exp_factor) / filter_value);
      position = position +1;

    end

    location = round((loops + 1) / 2);

    realv = median(values(1,:));
    imagv = median(values(2,:));

    oloop_output(hits(ii+1)+1) = realv + 1i*imagv;
    
  end
end

% UTILITY FUNCTIONS
function answer = floor_to_pow2(x)
  for ii = 1:x
      if ii<= x
          answer = ii;
      end
      ii = 2*ii;
  end
  answer = answer/2;
end

function answer = gcd(a,b)
  if (mod(a,b) == 0) 
    answer =  b;
  else 
    answer = gcd(b, mod(a,b)); 
  end
end

function answer = timesmod(x, a, n) 
 answer = round(mod((x * a),n));
end


function x = mod_inverse(a, m) 
	m0 = m; 
	y = 0;
	x = 1;
    a0 = a;
	if m == 1 
    x = 0;
  end

	while (a > 1)

    q = floor(double(a)/double(m)); 
    % disp("a/m");
    % disp(floor(a/m));
    % fprintf("a/m = %f\n",double(a)/double(m));
    % fprintf("q=%d, a=%d, m=%d\n",q,a,m);
    
		t = m ;
 
		m = mod(a,m);
		a = t ;
		t = y ;

		y = x - q * y ;
    x = t ;
    % fprintf("q=%d, a=%d, m=%d, x=%d, t=%d\n",q,a,m,x,t);
    
  end
  
	if x < 0
		x = x + m0 ;
  end
  % x = mod(x,m0);

end

function output = find_largest_indices(num,samples)

  % disp("log(samples)");
  % disp(log(samples));

  % disp("sum(log(samples)>3)");
  % disp(sum(log(samples)>3));

  % temp = sum(abs(samples)>2000);

  % if num > temp
  %   num = temp;
  % end

  if num > length(samples)
    num = length(samples);
  end
  

  [sorted,I] = sort(samples,'descend');
  % disp(size(I));
  % disp(num)
  % disp("sorted(1:num)");
  % disp(sorted(1:num));
  output = I(1:num)-1;
end 
