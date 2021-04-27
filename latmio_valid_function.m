function [Rlatt,CC,PL,SW,EG] = latmio_valid_function(R,ITERATIONS,D,pause)
%LATMIO_VALID_FUNCTION     Lattice with preserved degree distribution,   
%                          with metrics for validation
%
%   [Rlatt,CC,PL,SW,EG] = latmio_valid_function(R,ITERATIONS,D,pause);
%
%   This function "latticizes" an undirected network, while preserving the
%   degree distribution. The function does not preserve the strength
%   distribution in weighted networks. It also pulls out metrics of the null
%   network every n iterations, which can be used to validate the stability
%   and quality of the null model with the used iterations.
%
%   Input:      R,            undirected (binary/weighted) connection matrix
%               ITERATIONS,   rewiring parameter
%                             (each edge is rewired approximately ITER times)
%               D,            distance-to-diagonal matrix
%
%   Output:     Rlatt,        latticized network in original node ordering
%               CC/PL/SW/EG,  network metrics pulled out every nth ITER
%
%
%   Requires:   Various scripts from the BCT to calculate network metrics.
%
%
%   References: Maslov and Sneppen (2002) Science 296:910
%               Sporns and Zwi (2004) Neuroinformatics 2:145
%
%   2007-2012
%   Mika Rubinov, UNSW
%   Jonathan Power, WUSTL
%   Olaf Sporns, IU
%   Lance Burn, UC

%   Modification History:
%   Jun 2007: Original (Mika Rubinov)
%   Apr 2008: Edge c-d is flipped with 50% probability, allowing to explore
%             all potential rewirings (Jonathan Power)
%   Feb 2012: limit on number of attempts, distance-to-diagonal as input,
%             count number of successful rewirings (Olaf Sporns)
%   Feb 2012: permute node ordering on each run, to ensure lattices are
%             shuffled across mutliple runs (Olaf Sporns)
%   Mar 2021: Motified script to pull out network metrics every nth
%             iteration for null model validation (Lance Burn)
%             Original randmio_und.m function available in the BCT.


% real network
Kii=sum(R~=0,2);
cyc3ii=diag((R.^(1/3))^3);
Kii(cyc3ii==0)=inf;             %if no 3-cycles exist, make C=0 (via K=inf)
Cc=cyc3ii./(Kii.*(Kii-1));        %real clustering coefficient
C = mean(Cc);

clear Cc

CC = [];
PL = [];
SW = [];
EG = [];


n=size(R,1);

% randomly reorder matrix
ind_rp = randperm(n);
R = R(ind_rp,ind_rp);

% create 'distance to diagonal' matrix
if nargin<3 %if D is not specified by user
    D=zeros(n);
    u=[0 min([mod(1:n-1,n);mod(n-1:-1:1,n)])];
    for v=1:ceil(n/2)
        D(n-v+1,:)=u([v+1:n 1:v]);
        D(v,:)=D(n-v+1,n:-1:1);
    end
end
%end create

[i,j]=find(tril(R));
K=length(i);
ITER=K*ITERATIONS;

% maximal number of rewiring attempts per 'iter'
maxAttempts= round(n*K/(n*(n-1)/2));
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
            %lattice condition
            if (D(a,b)*R(a,b)+D(c,d)*R(c,d))>=(D(a,d)*R(a,b)+D(c,b)*R(c,d))
                R(a,d)=R(a,b); R(a,b)=0;
                R(d,a)=R(b,a); R(b,a)=0;
                R(c,b)=R(c,d); R(c,d)=0;
                R(b,c)=R(d,c); R(d,c)=0;
                
                j(e1) = d;          %reassign edge indices
                j(e2) = b;
                eff = eff+1;
                break;
            end %lattice condition
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
    
end %iterations

% lattice in node order used for latticization
Rrp = R;
% reverse random permutation of nodes
[~,ind_rp_reverse] = sort(ind_rp);
Rlatt = Rrp(ind_rp_reverse,ind_rp_reverse);
