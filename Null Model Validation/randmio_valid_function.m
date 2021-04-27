function [R,CC,PL,SW,EG]=randmio_valid_function(R, ITERATIONS, pause)
%RANDMIO_VALID_FUNCTION     Random graph with preserved degree
%                           distribution, with metrics for validation
%
%   R = randmio_valid_function(R,ITERATIONS,pause);
%   [R CC PL SW EG]=randmio_valid_function(R, ITERATIONS, pause);
%
%   This function randomizes an undirected network, while preserving the
%   degree distribution. The function does not preserve the strength
%   distribution in weighted networks. It also pulls out metrics of the null
%   network every n iterations, which can be used to validate the stability
%   and quality of the null model with the used iterations.
%
%   Input:      R,            undirected (binary/weighted) connection matrix
%               ITERATIONS,   rewiring parameter
%                             (each edge is rewired approximately ITER times)
%               pause,        metrics are recorded every nth iteration. 
%               
%
%   Output:     R,            randomized null model network
%               CC/PL/SW/EG,  network metrics pulled out every nth ITER
%
%
%   Requires:   Various scripts from the BCT to calculate network metrics.
%
%
%   References: Maslov and Sneppen (2002) Science 296:910
%
%
%   2007-2021
%   Mika Rubinov, UNSW
%   Jonathan Power, WUSTL
%   Olaf Sporns, IU
%   Lance Burn, UC

%   Modification History:
%   Jun 2007: Original (Mika Rubinov)
%   Apr 2008: Edge c-d is flipped with 50% probability, allowing to explore
%             all potential rewirings (Jonathan Power)
%   Mar 2012: Limit number of rewiring attempts, count number of successful
%             rewirings (Olaf Sporns)
%   Mar 2021: Motified script to pull out network metrics every nth
%             iteration for null model validation (Lance Burn)
%             Original randmio_und.m function available in the BCT.


Lnii = weight_conversion(R, 'lengths');
Dii = distance_wei(Lnii);
PL = charpath(Dii,0,0);       %real path length

CC = [];
PL = [];
SW = [];
EG = [];

n=size(R,1);
[i,j]=find(tril(R));
K=length(i);
ITER=K*ITERATIONS;

% maximal number of rewiring attempts per 'iter'
maxAttempts= round(n*K/(n*(n-1)));
% actual number of successful rewirings
eff = 0;

for iter=1:ITERATIONS
    att=0;
    while (att<=maxAttempts)                                     %while not rewired
        while 1
            e1=ceil(K*rand);
            e2=ceil(K*rand);
            while (e2==e1)
                e2=ceil(K*rand);
            end
            a=i(e1); b=j(e1);
            c=i(e2); d=j(e2);
            
            if all(a~=[c d]) && all(b~=[c d])
                break           %all four vertices must be different
            end
        end
        
        if rand>0.5
            i(e2)=d; j(e2)=c; 	%flip edge c-d with 50% probability
            c=i(e2); d=j(e2); 	%to explore all potential rewirings
        end
        
        %rewiring condition
        if ~(R(a,d) || R(c,b))
            R(a,d)=R(a,b); R(a,b)=0;
            R(d,a)=R(b,a); R(b,a)=0;
            R(c,b)=R(c,d); R(c,d)=0;
            R(b,c)=R(d,c); R(d,c)=0;
            
            j(e1) = d;          %reassign edge indices
            j(e2) = b;
            eff = eff+1;
            break;
        end %rewiring condition
        att=att+1;
    end %while not rewired
    
    
    if rem(iter,pause)==0
        %CC
        CCa = clustering_coef_wu(R);
        CCm = mean(CCa);
        CC = [CC CCm];
        
        %PL
        Lr = weight_conversion(R, 'lengths');
        Dr = distance_wei(Lr);
        PLrand = charpath(Dr,0,0);
        PL = [PL PLrand];
        
        %SW
        SWr=CCm/PLrand;
        SW = [SW SWr];
        
        %EG
        adjM_nrm = weight_conversion(R, 'normalize'); %rescale all edge weights to [0,1]
        EGccc = efficiency_wei(adjM_nrm); %normalization = binarization + ??
        EG = [EG EGccc];
    end
    
end


end %iterations
