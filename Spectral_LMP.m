function [y] = Spectral_LMP(Q,d,k,theta,x) 
y=x;
    for i = 1:k
    v = Q(:, i);                             
     alpha = 1 - theta/d(i);
    y = y - alpha * v * (v' * y);    
    end

end