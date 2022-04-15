function prob = SedumiToMosek(At,b,c,K)
    %Convert data in Sedumi Standard Primal form to Mosek 
    %Update: 04/15/2022
    %*******Important********
    %We only look consider the cone of free(K.f) and semidefinite(K.s)
    %The reason we seperate homogeneous case and unhomogeneous case is just
    %to precipitate the process
    
    %Consider 2 cases. 
    %1. K.s is homogeneous partition. 
    %2. K.s is not homogeneous.
    
    Homo = ChechHomo(K);
    NumOfPSD = length(K.s);
    NumOfFree = K.f;
    [hei,wei] = size(At);
    m = hei; %number of constraints
    At_PSD = At(:,NumOfFree+1:end);
    if Homo
        W_small = K.s(1); %dimension of each PSD
        %w = (K.s(1))^2;%Number of variables in each PSD
    end    
    
    
    
    %Objective
    c_f = c(1:NumOfFree); %free part of c
    c_s = c(NumOfFree+1:end); %PSD part of c
    prob.c = c_f;
    prob.blc = b;
    prob.buc = b;
    prob.blx = -inf*ones(1,NumOfFree);
    prob.bux = [];
    
    % Dimensions of PSD variables
    prob.bardim = K.s;    
    
    if Homo %fast computation
        [IndSym,IndDiag] = SymmetricIndices(W_small);
        dx = W_small*(W_small+1)/2;
        [row,col]=ind2sub([W_small,W_small],IndSym);
        
        c_rows=repmat(row,NumOfPSD,1);
        c_cols=repmat(col,NumOfPSD,1);
        c_vals=zeros((W_small+1)*W_small/2*NumOfPSD,1);
        c_IdxPsd=kron((1:NumOfPSD)',ones(dx,1)) ;
        start = 1; 
        InterestedIdx = repmat(IndSym,NumOfPSD,1);
        OffSet = kron(W_small^2*(0:NumOfPSD-1)',ones(dx,1));
        InterestedIdx = InterestedIdx+OffSet;
        %InterestedIdx = InterestedIdx + NumOfFree;
        %c_vals = c(InterestedIdx);
        c_vals = c_s(InterestedIdx);
        
        [r,c,v] = find(c_vals); 
        
        prob.barc.subj = c_IdxPsd(r);
        prob.barc.subk = c_rows(r);
        prob.barc.subl = c_cols(r);
        prob.barc.val = v; 
        
   
    else %Not Homegeneous we need to use loop
        c_rows = {};
        c_cols = {};
        c_vals = {};
        c_IdxPsd = {};
        OffSet = 0;
        for i = 1:NumOfPSD
            [IndSym,IndDiag] = SymmetricIndices(K.s(i));
            [row,col]=ind2sub([K.s(i),K.s(i)],IndSym);
            dx = K.s(i)*(K.s(i)+1)/2;
            c_rows{i} = row;
            c_cols{i} = col;
            c_IdxPsd{i}= i*ones(dx,1) ;
            start = 1; 
            InterestedIdx = IndSym;
            if i ~= 1
                OffSet = OffSet + K.s(i-1)^2;
            end
            InterestedIdx = InterestedIdx+OffSet;
            c_vals{i} = c_s(InterestedIdx);
        end
        
        IDXPSDs = vertcat(c_IdxPsd{:});
        ROWS = vertcat(c_rows{:});
        COLS = vertcat(c_cols{:});
        VALS = vertcat(c_vals{:});
        
        [r,c,v] = find(VALS);
        prob.barc.subj  = IDXPSDs(r);
        prob.barc.subk = ROWS(r);
        prob.barc.subl = COLS(r);
        prob.barc.val = v;
        
    end
    
    %constraint
    
    %linear part
    if NumOfFree ~= 0
        At_Free = At(:,1:NumOfFree);
        [r,c,v]=find(At_Free);
        prob.a = sparse(r,c,v,m,NumOfFree);
%         a_num = kron(1:m,ones(1,NumOfFree));
%         a_col = repmat(1:NumOfFree,1,m);
%         a_val = reshape(At(:,1:NumOfFree)',1,[]);
%         prob.a = sparse(a_num,a_col,a_val,m,NumOfFree);
    else
        prob.a = sparse([], [], [], m, 0); 
    end
    
    
    %PSD part
    if Homo  %fast computation
         A_rows=repmat(repmat(row,NumOfPSD,1),m,1);
         A_cols=repmat(repmat(col,NumOfPSD,1),m,1);
         A_vals=zeros(dx*NumOfPSD*m,1);

         A_IdxPsd = repmat(kron((1:NumOfPSD)',ones(dx,1)),m,1);
         A_InxCons= kron((1:m)',ones(dx*NumOfPSD,1));

         [r,c,v]=find(reshape(At_PSD(:,InterestedIdx)',[],1));    
         
         prob.bara.val = v;
         
         prob.bara.subi = A_InxCons(r);% Which constraint (i)
         prob.bara.subj = A_IdxPsd(r);
         prob.bara.subk = A_rows(r);
         prob.bara.subl = A_cols(r);         
                  
    else
        %We store the indices in the first constrait and use those
        %indices to get the value in other constraints
        A1_InxCon = {};
        A1_IdxPsd = {};
        A1_rows = {};
        A1_cols = {};
        InterestedIdxs = [];
        OffSet = 0;
        for i = 1:NumOfPSD
            [IndSym,IndDiag] = SymmetricIndices(K.s(i));
            [row,col]=ind2sub([K.s(i),K.s(i)],IndSym);
            dx = K.s(i)*(K.s(i)+1)/2;
            A1_rows{i} = row;
            A1_cols{i} = col;
            A1_IdxPsd{i} = i*ones(dx,1);
            A1_InxCon{i} = ones(dx,1);
            
            InterestedIdx = IndSym;
            if i>1
                OffSet = OffSet + K.s(i-1)^2;
            end
            InterestedIdx = InterestedIdx+OffSet;
            InterestedIdxs = [InterestedIdxs;InterestedIdx];
            
        end
        
        A1_r = vertcat(A1_rows{:});
        A1_c = vertcat(A1_cols{:});
        A1_IdPSDs = vertcat(A1_IdxPsd{:});
        A1_InxCons = vertcat(A1_InxCon{:});
        
        
        
        ROWS = repmat(A1_r,m,1);
        COLS = repmat(A1_c,m,1);
        IDXCONS = kron([1:m]',A1_InxCons);
        IDXPSDS = repmat(A1_IdPSDs,m,1);
        
        [r,c,v] = find(reshape(At_PSD(:,InterestedIdxs)',[],1));
        
        prob.bara.subk = ROWS(r);
        prob.bara.subl = COLS(r);
        prob.bara.subi = IDXCONS(r);
        prob.bara.subj = IDXPSDS(r);
        prob.bara.val = v;

        
    end
   
end