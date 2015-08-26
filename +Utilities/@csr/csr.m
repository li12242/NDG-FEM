classdef csr
    properties
        val     % val
        ser     % serial
    end %properties public
    properties(GetAccess=private)
        nRow    % row No.
        nCol    % col No.
        totalNum    % No. of value in the object
        localNum    % No. of value already stored in the object
    end %properties private
    
    methods
        function obj = csr(row, col, totalNum)
            obj.nRow = row;
            obj.nCol = col;
            obj.val = zeros(totalNum,1);
            obj.ser = zeros(totalNum,2);
            obj.localNum = 0;
            obj.totalNum = totalNum;
        end
        
        function obj = setVal(obj, irow, icol, val)
            assert(obj.localNum < obj.totalNum, ...
                'The csr matrix is full, no space left')
            obj.localNum = obj.localNum +1;
            obj.ser(obj.localNum, :) = [irow, icol];
            obj.val(obj.localNum) = val;
        end
        
%         function getVal(obj,irow,icol)
%         end
    end
end