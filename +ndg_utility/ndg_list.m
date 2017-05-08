classdef (HandleCompatible = true) ndg_list
    %NDG_LIST 生成二维链表数组
    %   二维链表按列存储
    %
    % Usages:
    %   生成包含5列向量的链表，其中每列元素个数分别为[2,2,5,4,5]
    %   list = ndg_list([2,2,5,4,5]);
    %   
    %   链表取值有两种方式，分别是按行列序号和按元素序号，如
    %   list(3)     % 返回第三个元素值，即第二列第一个元素;
    %   list(2, 4)  % 返回第二列第四个元素;
    %   注意，若是第4行只有1个元素，那么顺延至第5行取其第一个元素.
    %
    %   链表赋值有两种方式，分别是
    %   list = ndg_list([2,2,5]);
    %   list(3) = 1:7   % 从第二个元素开始，将1:7分别赋值到各个元素中;
    %   list(2,2) = 1:6 % 从第二行第2列开始，将1:6分别赋值到各个元素中;
    
    properties
        ndim; % 链表列数
        dim; % 各列元素数
        tol; % 元素总数
        var; % 元素值
    end
    
    methods
        function obj = ndg_list(elenum)
            obj.ndim = numel(elenum);
            obj.tol = 0;
            obj.dim = zeros(obj.ndim, 1);
            for n = 1:obj.ndim
                obj.dim(n) = elenum(n);
                obj.tol = obj.tol + elenum(n);
            end
            obj.var = zeros(obj.tol, 1);
        end
        
        function disp(obj)
            sk = 1;
            for n = 1:obj.ndim
                disp(['list(:, ', num2str(n), ')']);
                disp( obj.var( sk:(sk+obj.dim(n)-1) ) );
                sk = sk + obj.dim(n);
            end
        end% func
        
        function b = subsref(a, s)
            switch s.type
                case '()'
                    num = numel(s.subs);
                    if(num == 1) 
                        b = a.var(s.subs{1}); % for b = a(4) type
                    elseif( num == 2 )
                        if ( s.subs{1} == ':' ) % for b = a(:, 2) type
                            sk = 1;
                            len = a.dim(s.subs{2});
                            for n = 1:(s.subs{2}-1)
                                sk = sk + a.dim(n);
                            end
                            ind = sk:(sk+len-1);
                            b = a.var(ind);
                        else
                            sk = s.subs{1}; % for b = a(1, 2) type
                            for n = 1:(s.subs{2}-1)
                                sk = sk + a.dim(n);
                            end
                            b = a.var(sk);
                        end
                    end
                case '{}'
                case '.'
            end% switch
        end
        
        function a = subsasgn(a, s, b)
            % change vector from row to colume
            if isrow(b) 
                b = b'; 
            end
            
            switch s.type
                case '()'
                    len = numel(b); % length of value
                    num = numel(s.subs); % dimension of index
                    if(num == 1) 
                        ind = s.subs{1}:(s.subs{1}+len-1);
                        a.var(ind) = b;
                    elseif( num == 2)
                        if ( s.subs{1} == ':' ) % for a(:, 2) = b type
                            sk = 1;
                            len = a.dim(s.subs{2});
                            for n = 1:(s.subs{2}-1)
                                sk = sk + a.dim(n);
                            end
                            ind = sk:(sk+len-1);
                            a.var(ind) = b;
                        else
                            sk = s.subs{1}; % for a(2) = b type
                            for n = 1:(s.subs{2}-1)
                                sk = sk + a.dim(n);
                            end
                            ind = sk:(sk+len-1);
                            a.var(ind) = b;
                        end
                    end
                case '{}'
                case '.'
            end
        end
    end
    
end

