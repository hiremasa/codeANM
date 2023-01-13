clear;
clc;
%%%%%%%%%%%%%%%%%%%%
alpha = 0.05;
load Student;
data = Student;
x1 = round(data(:,31)/1);
indVal = [1,2,4,5,6,17,20,22]; % individual causes

%-------------check the two examples present in fig.5
k1 = data(:,1)*10 + data(:,2) + 1;  % x1 and x2
x1x2_p_val = [ANMs(k1, x1, alpha, 0),ANMs(x1, k1, alpha, 0)]
x1x2_d_val = MaxMI(k1,[data(:,1),data(:,2)],x1)
k2 = data(:,15)*10 + data(:,17) + 1; % x15 and x17
x15x17_p_val = [ANMs(k2, x1, alpha, 0),ANMs(x1, k2, alpha, 0)]
x15x17_d_val = MaxMI(k2,[data(:,15),data(:,17)],x1)

if x1x2_p_val(1)> alpha && x1x2_p_val(2) < alpha && x1x2_d_val<0.9
    combined_cause =  [1,2] % output combined_cause
end

if x15x17_p_val(1)> alpha && x15x17_p_val(2) < alpha && x15x17_d_val<0.9
    combined_cause =  [15,17] % output combined_cause
end

%-------------find all combined causes-start
%-------------which can be call several times
%-------------due to 'ANMs' is a random function with stochastic optimization
for i = 2
    M = nchoosek(1:30,i);
    [n,m]=size(M);
    for p = 1:n
        if  isempty(setdiff(M(p,:),intersect(M(p,:),indVal)))
            k = 0;
            for q = 1:m
                k = k + (data(:,M(p,q))+1)*10^(q-1);
            end
            resultM =  [ANMs(k, x1, alpha, 0),ANMs(x1, k, alpha, 0)];
            if resultM(1)< alpha && resultM(2) > alpha
                false_result_combined_cause = M(p,:)   % output any false results if exists
            end
            if resultM(1)> alpha && resultM(2) < alpha
                combined_cause =  M(p,:) % output combined_cause
            end
        end
    end
end
%-------------find all combined causes-end