function [ind, gSquare] = Conditional_independent_test(x, y, z, alpha)
% This is a function about independence test two groups variables.
% in this function, we use a group values same as in A but without
% repetition to replace A. It reduces the influence of singular value which
% far way from expection for A.
% tic
% rand('seed',1)
% z = ceil(rand(500,6) *10);
% x = z(:,3);
% z = z(:,[1:2,4:6]);
% y = x + ceil(rand(500,1) *10)-5;
% alpha = 0.05;

    len.x       = length(x);
    len.y       = length(y);
    [len.z, nz] = size(z);
    if len.x~= len.y || len.x~= len.z
        disp(len)
        error('the length among three input variables do not equal!');
    end
    
%%% find out no repetitions values in x and y. %%%
    
    uni.x       = unique(x);
    uni.y       = unique(y);    
    len.ux      = length(uni.x);
    len.uy      = length(uni.y);
    
    row.lz      = zeros(1,nz);
    row.mz      = zeros(2,nz);
    row.uz      = cell(1,nz);
    row.uz2     = zeros(len.z,nz);
    for i = 1:nz
        [row.uz{i}, ~, row.uz2(:,i)] = unique(z(:,i));
        row.mz(1,i)                  = max(row.uz{i});
        row.mz(2,i)                  = min(row.uz{i});
        row.lz(i)                    = length(row.uz{i});
    end
    len.uz      = prod(row.lz);
    C.uz        = zeros(len.uz, 1);
    clear i
    
%%% for matrix A. A_m*n*o(i,j,k) = 1  <=>  A_(i + (j-1)*m + (k-1)*(m*n)) = 1 
%%% <=> A_(p) = 1 that p = sum([i,j-1,k-1] .* [1, m, m*n]).  
%%% for the high dimension data set z, we makes it as a array newz. In newz
%%% , newz(i) = sum([z(1),z(2:end)-1] .*temp);
    temp         = cumprod(row.lz);
    temp         = [1, temp(1:end-1)];
    newz         = zeros(1,len.x);
    row.uz2      = row.uz2 - 1;
    row.uz2(:,1) = row.uz2(:,1) + 1;
    if ~all(row.uz2>0)
        row.uz2
        error('there is a element in z less than 0')
    end
    for i = 1:len.x
        newz(i)       = sum(row.uz2(i,:) .*temp);
        C.uz(newz(i)) = C.uz(newz(i)) + 1;
    end
    clear i temp
    
    uni.z       = unique(newz);
    C.uz        = C.uz(uni.z);
    len.uz      = length(C.uz);
    
    C.uxz       = zeros(len.ux, len.uz);
    C.uyz       = zeros(len.uy, len.uz);
    C.uxyz      = zeros(len.ux, len.uy, len.uz);
    
%%% Find out the distribution among x, y and z, and the their joint distribution. %%%
    row.z     = 1:len.z;
    row.ux    = 1:len.ux;
    row.uy    = 1:len.uy;
    row.uz    = 1:len.uz;
    
    num.zero  = 0;
    gSquare   = 0;
    for i = row.uz 
        temp       = row.z(newz == uni.z(i));
        set.x      = x(temp);
        set.y      = y(temp);
        
        for j = row.uy
            C.uyz(j,i)  = sum(set.y == uni.y(j));
        end
        clear j
        
        for j = row.ux
            tempx       = set.x == uni.x(j);
            C.uxz(j,i)  = sum(tempx);
            
            tempxz      = temp(tempx);
            if ~isempty(tempxz)
                tempsy      = y(tempxz);
                for k = row.uy
                    C.uxyz(j,k,i) = sum(tempsy == uni.y(k));
                    if C.uxyz(j,k,i) == 0
                        num.zero  = num.zero + 1;
                    else
                        tempxyz   = C.uxz(j,i) * C.uyz(k,i)/ C.uz(i);
                        gSquare   = gSquare + C.uxyz(j,k,i) * log(C.uxyz(j,k,i)/tempxyz);
                        clear tempxyz
                    end
                end
                clear k
            end
            clear tempx tempxz
        end
        clear j temp set        
    end
    
    temp      = prod(row.mz(1,:) - row.mz(2,:) +1);
    
    num.zero  = (max(x) - min(min(x),0)+1) * (max(y) - min(min(y),0)+1) * temp...
        - (len.ux * len.uy * len.uz) + num.zero;
    clear temp
    
    if num.zero < 0
        error('number of zero is less than 0.')
    end
    
    temp      = prod(row.mz(1,:) - row.mz(2,:));
    
    df        = (max(x) - min(x)) * (max(y) - min(y)) * temp - num.zero;
    df        = max(1, df);
    clear temp
    
    critical  = chi2inv(1 - alpha, df);    
    if 2*gSquare <= critical
        ind       = true;
    else
        ind       = false;
    end
    clear df critical len Max C num
%     toc
end