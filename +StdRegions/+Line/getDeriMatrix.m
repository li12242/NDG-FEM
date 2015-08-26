function Dr = getDeriMatrix(r)
% get Derivative Matrix
Dr = Polylib.Dglj(r);
Dr = Dr';
end