% Tim Truster
% Stabilized DG for elasticity problem, sectors are still portions of the
% element, since the Tau functions were already general enough.
% 05/16/15


%% Global declarations
% Place here any commands that should be executed
% whenever the element routine is called.
if ~exist('PSPS','var')
PSPS = 'n'; %plane stresS/straiN flag
end
nitvms = 1;2;
if nitvms == 1 %VMS
pencoeff = 4;1;
elseif nitvms == 2 %Nitsche
pencoeff = 1;4;2;
else %RFB
pencoeff = 1;4;2;
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

%%
    case {3,6,21} % Stiffness and internal force vector (REQUIRED)
        
        [ElemK,ElemF,hr] = L_Elem8_2d03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,MateT,isw,elem,numel,numSI,PSPS,lintt3,lintq4,lintt6,lintq9,nitvms,pencoeff);
        

%%
    case 51

        [ElemS,ElemE] = L_Elem8_2d51(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,PSPS,iprob,lintt3,lintq4,lintt6,lintq9);
        
        
%%
    case 53

        [ElemS,ElemE] = L_Elem8_2d53(mateprop,xl,ElemFlag,ndf,ndm,nst,nel,nen,PSPS,iprob,lintt3,lintq4,lintt6,lintq9);
        
        
end %Task Switch
