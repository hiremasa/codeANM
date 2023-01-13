function [ind, gSquare] = Independent_test(x, y, alpha)
% This is a function about independence test two groups variables.
% in this function, we use a group values same as in A but without
% repetition to replace A. It reduces the influence of singular value which
% far way from expection for A.

    len.x     = length(x);
    len.y     = length(y);
    if len.x ~= len.y
        disp(len);
        error('the length between x and y do not equal!')
    end
    
% find out no repetitions values in x and y
    
    uni.x     = unique(x);
    uni.y     = unique(y);
    
    len.ux    = length(uni.x);
    len.uy    = length(uni.y);
    C.ux      = zeros(len.ux, 1);
    C.uy      = zeros(len.uy, 1);
    C.uxy     = zeros(len.ux, len.uy);
    
%     find out the distribution of x and y, and the their joint distribution.
    row.x     = 1:len.x;
    row.ux    = 1:len.ux;
    row.uy    = 1:len.uy;
    
    for i = row.uy
       C.uy(i)    = sum(y == uni.y(i));
    end
    clear i
    
    num.zero  = 0;
    gSquare   = 0;
    for i = row.ux    
        temp       = row.x(x == uni.x(i));
        C.ux(i)    = length(temp); 
        
        set.y      = y(temp);
        for j = row.uy
            C.uxy(i,j)  = sum(set.y == uni.y(j));
            if C.uxy(i,j) == 0
                num.zero  = num.zero + 1;
            else
                temp       = C.ux(i) * C.uy(j)/ len.x;
                gSquare    = gSquare + C.uxy(i,j) * log(C.uxy(i,j)/ temp);
                clear temp
            end
        end
        clear j temp
    end
    
    num.zero  = (max(x) - min(min(x),1)+1) * (max(y) - min(min(y),1)+1) - (len.ux * len.uy) + num.zero;
    if num.zero < 0
        error('number of zero is less than 0.')
    end
    
    df        = (max(x) - min(x)) * (max(y) - min(y)) - num.zero;
    df        = max(1, df);
    critical  = chi2inv(1 - alpha, df);    
    if 2*gSquare <= critical
        ind       = true;
    else
        ind       = false;
    end
    gSquare = (2*gSquare)/critical;
    clear df critical len Max C num
end