%Test several instances of stable set problem

clc; clear;
addpath('Module\');
for pob = 2:8
    for NUM = 1:20
        filename = ['n30p',num2str(pob),'-',num2str(NUM)];
        Adj = readmatrix(['ErdosRenyiGraph\',filename,'.txt']);%adjacency matrix 
        n = width(Adj);
        [At_sdp,b_sdp,c_sdp,K_sdp]= FormStabelSetProblem_DualSDP(Adj);
        
        %At_sdp = full(At_sdp); b_sdp = full(b_sdp); c_sdp = full(c_sdp);
        At_sdp_copy = At_sdp;
        c_sdp_copy = c_sdp;
        OBJ_Outer_All = [];    
        for dx = [1,2,5] %partition

            J = ones(n,n);
            P = dx*ones(1,n/dx);
            opts.bfw = 1; opts.nop= n/dx; opts.dual = 0;

            [E,Comb,D]= ConstrOPR(P); %Construct Operator and the Combination
            NumOfComb = height(Comb);

            [Anew, bnew, cnew, Knew, info] = factorwidth_general(At_sdp,b_sdp,c_sdp,K_sdp,opts);
            Anew_copy = Anew;
            Acon = Anew(:,Knew.f+1:end);
            Acon_copy = Acon;
            U = eye(n);
            OBJ_Outer = [];
            Indices =BIGPSDposition(n,dx);
            IndicesAll =BIGPSDpositionAll(n,dx);
            len = (2*dx)*(2*dx+1)/2;
            for iter = 1:6
                %Mosek
                %prob1 = convert_sedumi2mosek(Anew,bnew,cnew,Knew);
                prob1 = SedumiToMosek(Anew,bnew,cnew,Knew);
                [rcode1, res1] = mosekopt('minimize info', prob1);
                OBJ_Outer = [OBJ_Outer,res1.sol.itr.pobjval];
                vecx = res1.sol.itr.barx;
                start = 1;
                X = zeros(n);
                tempx = zeros(n);
                start = 1;
                for num = 1:NumOfComb
                    X(Indices(num,:)) = X(Indices(num,:)) + vecx(start:start+len-1)';
                    start = start + len;
                    tempx(Indices(num,:)) = 0; %reset
                end  
                X = X+tril(X,-1)';
                X = U'*X*U;    
                U = chol(X);
                Acone=UpdateCone(Comb,NumOfComb,D,E,U);
                Anew(1:n^2,1+n^2+1:end) = Acone;
            end
            OBJ_Outer_All = [OBJ_Outer_All;OBJ_Outer];
        end
        
        %save(['SeveralRandomGraphs\\Result\\BSDD\\Outer\\',filename,'.mat'],'OBJ_Outer_All');
    end
end

%%
