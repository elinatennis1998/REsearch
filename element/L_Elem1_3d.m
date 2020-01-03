% Tim Truster
% 04/30/2013
%
% -Example 3D linear elastic solid element
% -Demonstrates the behavior of shifting dof ordering in FEAP
% -In isw=1 case, the element turns off all dofs that exceed ndm=3 because
%  they are not used for pure-displacement.
% -The same solution will be obtained, regardless of DG elements in the
%  mesh, if the ordering of dofs is changed in MatTypeTable:
%      MatTypeTable = [1; 8; 0]; 
%      MatTypeTable = [1; 8; 0; 3; 1; 2; 0]; 
%      MatTypeTable = [1; 8; 0; 3; 2; 1; 0]; 
%      MatTypeTable = [1; 8; 0; 1; 3; 2; 0];


%% Global declarations
% Place here any commands that should be executed
% whenever the element routine is called.


%% Task Switch - Perform various element functions
switch isw 
%%    
    case 1 % Setup up elemental degrees of freedom (OPTIONAL)
        
% Purpose: Rearrange order of dofs on element and/or flag certain
%          dofs to be presecribed away
% Called by: pmatin.m
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end

%%
    case 3 % Stiffness and internal force vector (REQUIRED)
        
        [ElemK,ElemF,hr] = L_Elem1_3d03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,ma,ieFEAP);
        

%%        
    case -1 % Boundary tractions (RECOMMENDED)
        
        ElemF = L_Elem1_3dm1(mateprop,nodeA,nodeB,elem,edge,traction,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen);
        
        
%%
    case 6 % Internal force vector (RECOMMENDED)
        
        [ElemK,ElemF,hr] = L_Elem1_3d03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,ma,ieFEAP);


%%
    case 21 % Stiffness matrix (RECOMMENDED)
        
        [ElemK,ElemF,hr] = L_Elem1_3d03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,ma,ieFEAP);


%%
    case 25 % Stress projection to nodes (RECOMMENDED)
        
        ElemS = L_Elem1_3d25(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nelS,nen,npstr,ma,ieFEAP);
        

%%        
    case 26 % Element Stress (OPTIONAL)
        
        ElemS = L_Elem1_3d26(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,nestr,ma,ieFEAP);
        

end %Task Switch
