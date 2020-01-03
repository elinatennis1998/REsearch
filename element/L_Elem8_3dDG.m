% Tim Truster
% 10/15/13
% DG for elasticity problem in 3D with conforming interfaces.
% The DG element input should simply be ix = [ixL ixR maDG] with ixL and
% ixR the connectivity from FormDG.

%% Global declarations
% Place here any commands that should be executed
% whenever the element routine is called.
nitvms = 1;
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
        
        if ndf > 3
            
            for i = 4:ndf
                lie(i,1) = 0;
            end
            
        end

%%
    case {3,6,21} % Stiffness and internal force vector (REQUIRED)
        
        [ElemK,ElemF,hr] = L_Elem8_3d03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,MateT,isw,elem,numel,numSI,nitvms,pencoeff);
        
        
end %Task Switch
