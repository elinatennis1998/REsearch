%Assembly od Internal and External Forces acring on Grain Edges
%Created 10/07/2019

if MatTypeTable(2,ma) == 10 
%     ForceI(:,ma) = ForceI(:,ma)+ElemFI_L; %+ ElemFI_R; %Interior Force Matrix for Left Elements
%     ForceE(:,ma) = ForceE(:,ma)+ElemFE_L; %+ ElemFE_R; %Exterior Force Matrix for Left Elements
%     ForceI(:,(ma+(nummat-nummatCG))) = ForceI(:,(ma+(nummat-nummatCG)))+ElemFI_R; %Interior Force Matrix for Right Elements
%     ForceE(:,(ma+(nummat-nummatCG))) = ForceE(:,(ma+(nummat-nummatCG)))+ElemFE_R; %Exterior Force Matrix for Right Elements
    MRDG_F_Ext(1:6,ma) = MRDG_F_Ext(1:6,ma) + ElemFE(1:6);
    MRDG_F_Int(1:6,ma) = MRDG_F_Int(1:6,ma) + ElemFI(1:6);
    MRDG_F_Ext(1:6,ma+nummat-nummatCG) = MRDG_F_Ext(1:6,ma+nummat-nummatCG) + ElemFE(7:12);
    MRDG_F_Int(1:6,ma+nummat-nummatCG) = MRDG_F_Int(1:6,ma+nummat-nummatCG) + ElemFI(7:12);

elseif MatTypeTable(2,ma) == 12
%     ForceI(:,ma) = ForceI(:,ma) + ElemFI_L; %Interior Force Matrix for Left Elements
%     ForceE(:,ma) = ForceE(:,ma) + ElemFE_L + ElemFE_R; %Exterior Force Matrix for Left Elements
    MRDG_F_Ext(1:6,ma) = MRDG_F_Ext(1:6,ma) + ElemFE;
    MRDG_F_Int(1:6,ma) = MRDG_F_Int(1:6,ma) + ElemFI;
elseif MatTypeTable(2,ma) == 11 %Solid Coarse Element Assembly
    MRDG_F_Ext(1:6,ma) = ElemFE; %Interior Force Matrix
    MRDG_F_Int(1:6,ma) = ElemFI; %Exterior Force Matrix
end

