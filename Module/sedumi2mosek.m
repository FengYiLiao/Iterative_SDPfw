function prob= sedumi2mosek(At,b,c,K)

   % n = sqrt(length(c));
    NumOfComb = length(K.s);
    [hei,wei] = size(At);
    m = hei;
    W_small = K.s(1);
    w = (K.s(1))^2;
    %[IndSym] = SymmetricIndices(W_small);
    % The scalar part, as in linear optimization examples
    prob.c = [];
    prob.a = sparse([], [], [], m, 0);          % 2 constraints, no scalar variables
    prob.blc = [b];                        % Bounds
    prob.buc = [b];
    
    
    % Dimensions of PSD variables
    prob.bardim = K.s;
    
    [IndSym,IndDiag] = SymmetricIndices(W_small);
    dx = W_small*(W_small+1)/2;
    %dx = sqrt(D(num))*(sqrt(D(num))+1)/2;
    [row,col]=ind2sub([sqrt(w),sqrt(w)],IndSym);
    c_rows=repmat(row,NumOfComb,1);
    c_cols=repmat(col,NumOfComb,1);
    c_vals=zeros((sqrt(w)+1)*sqrt(w)/2*NumOfComb,1);
    c_IdxPsd=kron((1:NumOfComb)',ones(dx,1)) ;
    start = 1; 
    InterestedIdx = repmat(IndSym,NumOfComb,1);
    OffSet = kron(W_small^2*(0:NumOfComb-1)',ones(dx,1));
    InterestedIdx = InterestedIdx+OffSet;
    c_vals = c(InterestedIdx);
%     for num = 1:NumOfComb
%         fprintf('%d\n',num);
%         ci = kron(EE(:,:,num)',EE(:,:,num)')*c;
%         ci_s = ci(IndSym);
%         c_vals(start:start+dx-1) = ci_s; 
%         %c_rows = []
%         %c_new(start:start+D(num)-1) = kron(EE(:,:,num)',EE(:,:,num)')*c;
%         %c_cmpt(:,num) =kron(EE(:,:,num)',EE(:,:,num)')*c;
%         %start = start + (D(num)+1)D(num);
%         start = start + dx;
%     end 
    
    prob.barc.subj = c_IdxPsd;
    prob.barc.subk = c_rows;
    prob.barc.subl = c_cols;
    prob.barc.val = c_vals;
    
    
    %constraint

    %A_rows= repmat(repmat())
    %A_rows=repmat(repmat(repmat(row,NumOfComb,1),m,1),,);
   % dx = (sqrt(w)+1)*sqrt(w)/2;
    
%     A_rows = [];
%     A_cols = [];
%     A_vals = [];
%     A_IdxPsd = [];
%     A_InxCons = [];
    
    A_rows=repmat(repmat(row,NumOfComb,1),m,1);
    A_cols=repmat(repmat(col,NumOfComb,1),m,1);
    A_vals=zeros(dx*NumOfComb*m,1);
    %A_IdxPsd= kron((1:m)',ones(dx*NumOfComb,1));
    A_IdxPsd = repmat(kron((1:NumOfComb)',ones(dx,1)),m,1);
    A_InxCons= kron((1:m)',ones(dx*NumOfComb,1));
    
%     InterestedAt = At(:,InterestedIdx);
%     InterestedAt = InterestedAt';
%    A_vals = reshape(At(:,InterestedIdx)',[],1);
    
     prob.bara.subi = A_InxCons;% Which constraint (i)
     prob.bara.subj = A_IdxPsd;
     prob.bara.subk = A_rows;
     prob.bara.subl = A_cols;
     prob.bara.val = reshape(At(:,InterestedIdx)',[],1);    
    
    
end