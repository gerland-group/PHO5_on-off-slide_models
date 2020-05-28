function ps = calc_stat_distr_state_reduction(W_gain)  
% Calculates the steady state distribution of the rate matrix W_gain
% Implementation of the state reduction algorithm aka GTH (Grassmann, Taksar, Heyman) algorithm:
% Works for an irreducible, recurrent Markov chain. W' (calculated below from W_gain) is the generator matrix of a continous-time Markov chain.
% The state reduction algorithm is numerically very stable, since it does not use subtractions.
% Use the matlab coder to generate C code (mex function) for >2x speedup.
    
    if any(any(W_gain-diag(diag(W_gain))<0))
        W_gain
        error('W_gain has non-diagonal entries < 0!')
    end
    
    W = W_gain - diag(sum(W_gain));  % if W_gain is already a W-type matrix, this doesn't change anything
    Q = W';
   
    q = zeros(1,size(W,1));  % unnormalized steady state
    
    n_max = size(Q,1);
    for n = n_max:-1:2  % reduce number of states
        n_reduced = n-1;
        temp_sum = sum(Q(n,1:n_reduced));
        Q(1:n_reduced,n) = Q(1:n_reduced,n) / temp_sum;  % calculate reduced generator / transition probabilities
        for j=1:n_reduced
            Q(1:n_reduced,j) = Q(1:n_reduced,j) + Q(1:n_reduced,n) * Q(n,j);
        end
    end
    
    q(1) = 1;  % with one state left the steady state vector is just 1
    
    for k = 2:n_max  % rebuild steady state vector
        k_reduced = k-1;
        q(k) = sum(q(1:k_reduced).*(Q(1:k_reduced,k))');
    end
    
    ps = q'/(sum(q));  % normalize
    
end
