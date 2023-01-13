clc;
clear;
alpha = 0.05;
load InflamData
data = InflamData;
x1 = data(:,7); 
x2 = data(:,8); 
y1 = data(:,1); 
y2 = data(:,2); 
y3 = data(:,3); 
y4 = data(:,4); 
y5 = data(:,5); 
y6 = data(:,6); 
fprintf('1) Use CCI to find combined cause \n\n');
for i=1:3
    M = nchoosek(1:6,i);
    [n,m]=size(M);
    for p = 1:n
        k = 0;
        for q = 1:m
            k = k + (data(:,M(p,q))+1)*10^(q-1);
        end 
        resultM =  [ANMs(k, x1, alpha, 0),ANMs(x1,k, alpha, 0)];
        if resultM(1)> alpha && resultM(2) < alpha && ~isempty(intersect(M(p,:),6)) && CMaxMI(k,data(:,M(p,:)),x1)<0.9
           combined_cause = M(p,:)
           p_val = resultM
           d_val = CMaxMI(k,data(:,M(p,:)),x1)
        end
    end
end
fprintf('2) Use MHPC to find combined cause \n\n');
combined_cause_detected_By_MMHPC = MHPC([],[])