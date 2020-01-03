% Two dimenisional plane stress/strain element

%% Global declarations
% Place here any commands that should be executed
% whenever the element routine is called.
if ~exist('PSPS','var')
PSPS = 'n'; %plane stresS/straiN flag
end
if ~exist('iprob','var')
iprob = 0;
end


%% Task Switch - Perform various element functions
switch isw 
%%    
    case 1 % Setup up elemental degrees of freedom (OPTIONAL)
        
% Purpose: Rearrange order of dofs on element and/or flag certain
%          dofs to be presecribed away
% Called by: pmatin.m
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
        istv = 7; % number of stresses per node for post-processing
        iste = 7; % number of stresses per element

%%
    case 3 % Stiffness and internal force vector (REQUIRED)
        
        [ElemK,ElemF,hr] = L_Elem1_2d03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,lintt3,lintq4,lintt6,lintq9,PSPS,iprob);
        

%%        
    case -1 % Boundary tractions (RECOMMENDED)
        
        ElemF = L_Elem1_2dm1(mateprop,nodeA,nodeB,elem,edge,traction,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen);
        
        
%%
    case 6 % Internal force vector (RECOMMENDED)
        
        [ElemK,ElemF,hr] = L_Elem1_2d03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,lintt3,lintq4,lintt6,lintq9,PSPS,iprob);


%%
    case 21 % Stiffness matrix (RECOMMENDED)
        
        [ElemK,ElemF,hr] = L_Elem1_2d03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,lintt3,lintq4,lintt6,lintq9,PSPS,iprob);


%%
    case 25 % Stress projection to nodes (RECOMMENDED)
        
        ElemS = L_Elem1_2d25(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelS,nen,npstr,PSPS,iprob);
        

%%        
    case 26 % Element Stress (OPTIONAL)
        
        ElemS = L_Elem1_2d26(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,nestr,PSPS,iprob);
        

%%
    case 50

        [ElemS,ElemE,ElemD,ElemV,ElemW] = L_Elem1_2d50(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,PSPS,iprob,lintt3,lintq4,lintt6,lintq9);
        
        
%%
    case 52

        [ElemS,ElemE,ElemV] = L_Elem1_2d52(mateprop,xl,ElemFlag,ndf,ndm,nst,nel,nen,PSPS,iprob,lintt3,lintq4,lintt6,lintq9);
        
        
end %Task Switch
