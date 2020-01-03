% Tim Truster
% 10/20/2015
% 3D CZM element, elastic stiffness only


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
    case {3,6} % Stiffness and internal force vector (REQUIRED)
        
        [ElemK,ElemF,hr] = L_Elem7_3d03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen);
        
        
end %Task Switch
