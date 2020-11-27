WITH_COMB  = false;
ALGORITHM1 = true;
VERBOSE    = false;
% TIMING     = false;

function oloop_output = outer_loop(origx,n, filter_timedo,filter_sizet,filter_freqdo, 
  filter_est_timedo,filter_est_sizet,filter_est_freqdo, B2,
  	  num,B,W_Comb,Comb_loops,loop_threshold,location_loops,
	  loops)

    permute_vals = zeros(1,loops);
    permuteb_vals = zeros(1,loops);

% complex_t *x_samp[loops];
% for(int i = 0; i < loops; i++){
%   if (i < location_loops)
%     x_samp[i] = (complex_t*)calloc(B, sizeof(*x_samp[i]));
%   else
%     x_samp[i] = (complex_t*)calloc(B2, sizeof(*x_samp[i]));
% }    
    x_samp_B = zeros(B,loops);
    x_samp_B2 = zeros(B2,loops);

  hits_found =0;

  if WITH_COMB
    pages_hit = num * (n/B) * num * (1 / W_Comb) * location_loops;
  else
    pages_hit = num * (n/B) * location_loops;
  end

 score = zeros(1,n);
 hits = zeros(1,n);
  
  PF_LOC = 0;
  G_LOC = 0;

%  BEGIN INNER LOOPS

  for ii=0:loops-1
    a = 0;
    b = 0;

    while gcd(a, n) ~= 1
      a = int(mod(randi(n),n));
    end

    ai = mod_inverse(a, n);

    permute_vals(ii+1) = ai;
    permuteb_vals(ii+1) = b;

    perform_location = int64(i < location_loops);

    if perform_location
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Filter cur_filter = filter;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cur_B =  B;
    else
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Filter cur_filter = filter_Est;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      cur_B =  B2;
    end

    J = zeros(1,num);

    [x_samp(ii+1),J] = inner_loop_locate(origx, n, cur_filter,
                      num, cur_B,
                      a, ai, b, J);

    if perform_location 
      [score, hits, hits_found] = inner_loop_filter_regular(J, n, num, cur_B,
                                  a, ai, b, loop_threshold,
                                  score, hits, hits_found);
    end

  
  end
  
  %END INNER LOOPS
  oloop_output = estimate_values(hits,hits_found,
      x_samp, loops,n, permute_vals, B, B2,
      filter_timedo,filter_sizet,filter_freqdo, 
      filter_est_timedo,filter_est_sizet,filter_est_freqdo,location_loops)
end


function [x_samp,J] = inner_loop_locate(origx,n, filter_timedo,filter_sizet,filter_freqdo,
                     num,B,a,ai,b){

  if mod(n,B) ~= 0
    disp("Warning: n is not divisible by B, which algorithm expects.\n");
  end
  
  x_sampt = zeros(1,n);

  index=b;
  for ii = 0:filter_sizet
    x_sampt[mod(ii,B)+1] += origx[index+1] * filter_timedo[ii+1];
    index = mod((index+ai),n);
  end

  x_samp = fft_recur(x_sampt);

  samples = abs(x_samp[:B]);

  J= find_largest_indices(num,samples);
}


function [score, hits, hits_found] = inner_loop_filter_regular(J,n,num,B,a,ai,b,loop_threshold,
                              score, hits, hits_found)

  for ii=0:num-1
    low = mod((int64(ceil((J(i) - 0.5) * n / B)) + n),n);
    high = mod((int64(ceil((J(i) + 0.5) * n / B)) + n),n);
    loc = timesmod(low, a, n);
    
    jj= low;
    while jj ~= high
      score(loc+1) = score(loc+1)+1;

        if score(loc)==loop_threshold 
          hits(hits_found+1)=loc;
          hits_found +=1;
          loc = mod((loc + a),n);
        end

      jj = mod((jj+1),n);
    end
  end

end

function oloop_output = estimate_values(hits,hits_found,
                x_samp, loops,n, permute_vals, B, B2,
                filter_timedo,filter_sizet,filter_freqdo, 
                filter_est_timedo,filter_est_sizet,filter_est_freqdo,location_loops)
  
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
      permuted_index= timesmod(permute_vals(jj), hits(ii),  n);
      hashed_to = (permuted_index / (n / cur_B));
      dista = mod(permuted_index , (n / cur_B));
     
      if dista > ((n/cur_B)/2) 
        hashed_to = mod((hashed_to + 1),cur_B);
        dista = (dista -(n/cur_B));
      end
      dista = mod((n - dista),n);

      filter_value = cur_filter_freqdo(dista);
      
      values(1)(position) = real(x_samp(jj)(hashed_to) / filter_value);
      values(2)(position) = imag(x_samp(jj)(hashed_to) / filter_value);
      position = position +1;

    end

    location = int64((loops + 1) / 2);

    values(1) = median(values(1));
    values(2) = median(values(2));
    
    realv = values(1)(location);
    imagv = values(2)(location);

    oloop_output(hits(i)) = realv + 1i*imagv;
    
  end
end