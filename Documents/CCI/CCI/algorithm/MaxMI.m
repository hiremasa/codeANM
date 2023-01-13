
function  dval = MaxMI(CC,data,y)
[~,m] = size(data);
plusVal = 0;
MaxG = 0;
for i=1:m-1
    M = nchoosek(1:m,i);
    [n2,m2]=size(M);
    for p = 1:n2
        k = 0;
        for q = 1:m2
            k = k + (data(:,M(p,q))+1)*100^(q-1);
        end
           gSquare = CIxy_d(k,y)+plusVal;
        if gSquare >= MaxG
            MaxG = gSquare;
        end
    end
end
dval = MaxG/CIxy_d(CC,y);
end
