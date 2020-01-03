% Boundary condition weak enforcement 
% Tim Truster 
% Created 9/7/2019 
% Last Modified 9/8/2019


%% Global declarations
% Place here any commands that should be executed
% whenever the element routine is called.
if ~exist('PSPS','var')
PSPS = 'n'; %plane stresS/straiN flag
end
if ~exist('DispFine','var')
    DispFine = 0; %plane stresS/straiN flag
elseif ~exist('DispCoarse','var')
    DispCoarse = 0;
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
        iste = 20; % number of stresses per element

%%
    case {3,6,21} % Stiffness and internal force vector (REQUIRED)
        
               
        if exist('IBP','var') && IBP == 1
            [ElemK,ElemF,hr] = L_Elem12_2d03BConly(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,...
                ndf,ndm,nst,nel,nen,MateT,isw,elem,numel,numSI,PSPS,lintt3,lintq4,lintt6,lintq9,...
                nitvms,pencoeff,BoundXmid,nummatCG,ma,MaterTypeNum,TFS,factorT,factorS,eRVE,BoundEps,BoundSig,...
                RegionsOnInterface,GrainXmid,DmatSachs,BoundSachs,DispFine);
        else
                [ElemK,ElemF,hr] = L_Elem12_2d03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,...
            ndf,ndm,nst,nel,nen,MateT,isw,elem,numel,numSI,PSPS,lintt3,lintq4,lintt6,lintq9,...
            nitvms,pencoeff,BoundXmid,nummatCG,ma,MaterTypeNum,TFS,factorT,factorS,eRVE,BoundEps,BoundSig,...
            RegionsOnInterface,GrainXmid,DmatSachs,BoundSachs,DispFine,ElemE,ElemS,bCrys);
        end


%%
    case 51

        [ElemS,ElemE] = L_Elem12_2d51(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,PSPS,iprob,lintt3,lintq4,lintt6,lintq9);
        
        
%%
    case 53

        [ElemS,ElemE] = L_Elem12_2d53(mateprop,xl,ElemFlag,ndf,ndm,nst,nel,nen,PSPS,iprob,lintt3,lintq4,lintt6,lintq9,MateT,isw,elem,numel,numSI,nitvms,pencoeff);

        
%%
    case 59 % MRDG forces
        
        if exist('MRDGF','var') && MRDGF == 2
            [ElemFEA,ElemFEB,ElemFEC,ElemFED,ElemFEE,ElemFEF,ElemFIA,ElemFIB,ElemFIC,hr] = ...
                L_Elem12_2d59(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,...
                ndf,ndm,nst,nel,nen,MateT,isw,elem,numel,numSI,PSPS,lintt3,lintq4,lintt6,lintq9,...
                nitvms,pencoeff,BoundXmid,nummatCG,ma,MaterTypeNum,TFS,factorT,factorS,eRVE,BoundEps,BoundSig,...
                RegionsOnInterface,GrainXmid,DmatSachs,BoundSachs,DispFine,IBP,bCrys);
            MRDG_F_ExtA(1:6,ma) = MRDG_F_ExtA(1:6,ma) + ElemFEA(1:6);
            MRDG_F_ExtB(1:6,ma) = MRDG_F_ExtB(1:6,ma) + ElemFEB(1:6);
            MRDG_F_ExtC(1:6,ma) = MRDG_F_ExtC(1:6,ma) + ElemFEC(1:6);
            MRDG_F_ExtD(1:6,ma) = MRDG_F_ExtD(1:6,ma) + ElemFED(1:6);
            MRDG_F_ExtE(1:6,ma) = MRDG_F_ExtE(1:6,ma) + ElemFEE(1:6);
            MRDG_F_ExtF(1:6,ma) = MRDG_F_ExtF(1:6,ma) + ElemFEF(1:6);
            MRDG_F_IntA(1:6,ma) = MRDG_F_IntA(1:6,ma) + ElemFIA(1:6);
            MRDG_F_IntB(1:6,ma) = MRDG_F_IntB(1:6,ma) + ElemFIB(1:6);
            MRDG_F_IntC(1:6,ma) = MRDG_F_IntC(1:6,ma) + ElemFIC(1:6);
        end
        
end %Task Switch
