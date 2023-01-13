function [skeleton, data] = Data_generator_double(times,skeleton1,samN)

%     Gnum.var      : the number of node.
%     Gnum.nootnode  : the number of root node.
%     2  : the number of in-degree of each child node.
%     skeleton(i,j)  : if the node j is the parent of node i. 1 means yes,
%     and 0 means not.
%     temp           : the temporary variables whice are only used in this
%     funcion.
%     data           : the final discrete data.
%     temp.value     : the noise.
%     Gnum.sam     : the number of samples.


%     skeleton     = Skeleton_generator(times);
     skeleton = skeleton1';
    data         = Data_genetor(times, skeleton,samN);
%     disp(['---------------------- All done! ',num2str(times),' times ----------------------'])
end

function skeleton = Skeleton_generator(times)

%     Gnum.var      : the number of variable.
%     Gnum.root  : the number of root variable.
%     2  : the number of in-degree of each child variable.
%     skeleton(i,j)  : if the variable j is the parent of variable i. 1 means yes,
%     and 0 means not.

    global Gnum
    
    if Gnum.root < 2
        error('The number of root variable must bigger than in-degree!');
    end
    
    skeleton       = zeros(Gnum.var);
    temp.row       = 1:Gnum.var;
    
    temp.inrow     = ones(1, Gnum.var);
    temp.ninrow    = floor((Gnum.var -Gnum.root) *(1-Gnum.single));
    temp.inrow2    = randperm(Gnum.var -Gnum.root) +Gnum.root;
    temp.inrow(temp.inrow2(1:temp.ninrow)) = 2;
    temp.inrow(1:Gnum.root)            = 0;
    temp.parent                            = cell(Gnum.var, 1);
    
    Gnum.inrow     = temp.inrow;
    
    for i = Gnum.root +1 : Gnum.var
        temp.pnum      = 0;
        while temp.pnum < Gnum.inrow(i)
            parent         = temp.row(randi(i -1));
            if skeleton(temp.row(i), parent) == 0
                prow = [];
                for j = 1:length(temp.parent{i})
                    temp.prow   = prow;
                    prow        = [temp.prow,temp.parent{temp.parent{i}(j)}];
                end
                if ~any(parent == temp.parent{parent}) && ~any(parent == prow)
                    skeleton(temp.row(i), parent) = 1;
                    temp.pnum           = temp.pnum +1;
                    temp.parent{i}      = [temp.parent{i}, parent];
                end
            end
        end
    end
    
    route        = ['Variable_',num2str(Gnum.var),'_Sample_',num2str(Gnum.sam),...
        '_Runs_',num2str(times)];
    save(['Data\Actual_skeleton\', route,'.dat'],'skeleton','-ASCII');
    
%     disp('Skeleton is generated!')
    clear i temp route array Gnum
    
end

function data = Data_genetor(times, skeleton,samN)
    global Gnum
    Gnum.alpha           = 0.05;
    Gnum.times           = 5;
    Gnum.var             = 15;
    Gnum.sam             = samN;
    Gnum.root            = 3;
    Gnum.single          = 0;
    [num.m, num.n] = size(skeleton);
    if num.m ~= num.n
        error('The skeleton is not square!\n')
    end
    
    temp.ninrow    = sum(skeleton, 2);
    temp.ninrow    = temp.ninrow(temp.ninrow >0);
    temp.root      = zeros(1,num.n);
    temp.row       = 1:num.n;
    for i = 1:num.n
        if all(skeleton(i, :) == 0)
            temp.root(i)   = 1;
        end
    end
    temp.child     = temp.row(temp.root == 0);
    temp.root      = temp.row(temp.root == 1);
    len.troot      = length(temp.root);
    len.tchild     = length(temp.child);
%     if len.troot ~= Gnum.root
%         error('The number of root node is error!\n');
%     end
    
    %%% Create root nodes
    num.n;
    Gnum.sam;
    data           = zeros(Gnum.sam, num.n);
    
    for i = 1:len.troot
        p.x                   = 0.8 * rand + 0.1;
        data(:, temp.root(i)) = binornd(2, p.x, Gnum.sam, 1);
        clear p
    end
    
    
    %%% Create children nodes
    data;
    for i = 1:len.tchild
        parent         = data(:, skeleton(temp.child(i),:) == 1);
        if size(parent,2) ~= temp.ninrow(i)
            error('The skeleton data is error!\n');
        end
        
        n           = size(parent, 2);
        if n>1
            state       = zeros(1,n);
            for j = 1:n
                state(j)    = length(unique(parent(:,j)));
            end
            p.CPT       = 0.4 + rand*0.2;
            CPT         = binornd(1, p.CPT, prod(state), 1);
            clear n state
            
            c.parent    = coding(parent);
            Y           = CPT(c.parent);
        else
            Y           = mod(parent, 2);
        end
        
                 N           = binornd(1, 0.2 * rand + 0.4, Gnum.sam, 1);
%         N = 0;
        Y           = Y+N;
        
        data(:, temp.child(i)) = Y;
        
        clear Y MPV parent len_N N
    end
    
    route        = ['Variable_',num2str(Gnum.var),'_Sample_',num2str(Gnum.sam),...
        '_Runs_',num2str(times)];
%    save(['Data\Discrete_data\',route,'.dat'],'data','-ASCII');
    
%     disp('Data set is generated!')
    clear temp i route
end



function c_p = coding(p)
    %%%% matrix/vector is coded with unique variable of itself. %%%%
    [m, n]      = size(p);
    state.len   = zeros(1, n);
    state.var   = cell(1, n);
    for i = 1:n
        state.var{i} = unique(p(:,i));
        state.len(i) = length(state.var{i});
        pos          = cell(1,n);
        for j = 1:state.len(i)
            pos{j}        = (p(:,i) == state.var{i}(j));
        end
        clear j
        for j = 1:state.len(i)
            p(pos{j}, i)  = j;
        end
        clear j pos
    end
    clear i
    
    c_p         = zeros(m, 1);
    temp        = n:-1:1;
    temp        = state.len(temp);
    temp        = cumprod(temp);
    temp        = [1, temp(1 : end-1)];
    temp        = sort(temp, 'descend');
    row         = ones(1, n);
    row(end)    = 0;
    for i = 1:m
        c_p(i)       = sum((p(i, :) - row) .* temp);
    end
    clear i temp row m n p state    
end