function p_val = ANMs(X, Y, level, type)
    
    Max_times         = 10;
    num.MPV           = min( max(Y)-min(Y), 20);
    
    [val.X, ~, row.X] = unique(X);
    val.Y             = (min(Y):max(Y))';
    
    len.X             = length(val.X);
    len.Y             = length(val.Y);
    
    if size(val.X,1) == 1 || size(val.Y,1) == 1
        MPV           = ones(length(val.X), 1)* val.Y(1);
        p_val         = 1;
    else
        p             = hist3([X,Y], {val.X, val.Y});
        
        MPV           = zeros(len.X, 1);
        temp.MPV      = cell(1,len.X);
        for i = 1 : len.X
            [~, row.p]         = sort(p(i,:), 'ascend');
            
            for j = 1 : len.Y
                if j ~= row.p(end)
                    p(i,j)     = p(i,j) + 1/(2* abs(j - row.p(end)));
                else
                    p(i,j)     = p(i,j) + 1;
                end
            end
            
            [~, row.p]         = sort(p(i,:), 'ascend');
            temp.MPV{i}        = row.p;
            MPV(i)             = val.Y(row.p(end));
            clear row.p i
        end
        
        row.Y         = MPV(row.X);
        val.residual  = Calculate_Residual(Y, row.Y, type);        
        len.residual  = length(unique(val.residual));
        
        if len.residual == 1
            display('Warning!! there is a deterministic relation between X and Y');
            p_val              = 1;
        else
            p_val              = Chi_square(val.residual, X, len.residual, len.X);
        end
        
        times         = 0;
        while (p_val < level) && (times < Max_times)
            for i = randperm(len.X)
                row.MPV        = cell(1, num.MPV + 1);
                p_val_1        = zeros(1, num.MPV + 1);
                p_val_2        = zeros(1, num.MPV + 1);
                for j = 1 : num.MPV + 1
                    row.MPV{j}                = MPV;
                    row.MPV{j}(i)             = val.Y(temp.MPV{i}(end - (j-1)));
                    row.Y                     = row.MPV{j}(row.X);
                    val.residual              = Calculate_Residual(Y, row.Y, type);
                    [p_val_1(j), p_val_2(j)]  = Chi_square(val.residual, X, len.residual, len.X);
                end
                
                [m, n]             = max(p_val_1);
                if m < 1e-3
                    [~, n]             = min(p_val_2);
                end
                
                MPV                = row.MPV{n};
                row.Y              = MPV(row.X);
                val.residual       = Calculate_Residual(Y, row.Y, type);
                p_val              = Chi_square(val.residual, X, len.residual, len.X);
                
                clear p_val_1 p_val_2 row.MPV j m n
            end
            times              = times + 1;
            clear i
        end
        if type == 0
            MPV                = MPV + round(mean(val.residual));
        end
%         [~, ~, row.MPV]           = unique(MPV);
%         MPV                       = MPV(row.MPV);
        clear p temp.MPV X Y row len val times Max_times num level
    end
end

function residual = Calculate_Residual(Y, rowY, type)
    if type == 0
        residual   = Y - rowY;
    else
        residual   = mod(Y - rowY, max(Y) - min(Y) + 1);
    end
end