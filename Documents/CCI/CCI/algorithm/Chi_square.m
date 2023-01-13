function [p_val, T]  = Chi_square(X, Y, lenX, lenY)
    [~, ~, X]     = unique(X);
    [~, ~, Y]     = unique(Y);
    X             = X - min(X);
    Y             = Y - min(Y);
    len.X         = lenX;
    len.Y         = lenY;
    row.X         = 0 : len.X - 1;
    row.Y         = 0 : len.Y - 1;
    len.newX      = length(X);
    
    if len.X == 1 || len.Y == 1
        p_val       = 1;
        T           = 0;
    else
        p           = hist3([X, Y], {row.X, row.Y});
        mp          = sum(p, 1);
        np          = sum(p, 2);
        null_mp     = sum(mp == 0);
        null_np     = sum(np == 0);
        
        matrix.one  = zeros(len.X, len.Y);
        matrix.two  = zeros(len.X, len.Y);
        for i = 1 : len.X
            for j = 1 : len.Y
                matrix.one(i,j)    = (np(i)* mp(j))/ len.newX;
                if matrix.one(i,j) > 0
                    matrix.two(i,j) = (p(i,j) - matrix.one(i,j))^2/ matrix.one(i,j);
                end
            end
        end
        T           = sum(sum(matrix.two));
        p_val       = 1 - chi2cdf(T, (len.X -1 -null_np)* (len.Y -1 -null_mp));
        clear i j matrix len row null_mp null_np np mp p X Y
    end
end