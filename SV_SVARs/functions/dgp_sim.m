function [y] = dgp_sim(n,T,p, B0, a, A, D)
% Creates a y process from given parameters, according to the Equations
% provided in the text

 %   k = n*p+1;
    y = zeros(T,n);

    y_init = repelem(a,1, p);

    for t=1:T
       
        epsilon = mvnrnd(zeros(n,1), squeeze(D(:,:,t)));
        
        ut = B0 \ epsilon';

        if t <= p
            y_t = y_init(:, t);
        else
            y_t = a;
            for j = 1:p
                y_t = y_t + A(:, :, j) * y(t - j,:)';
            end     

        end

        y_t = y_t + ut;

        y(t, :) = y_t;

    end 



end