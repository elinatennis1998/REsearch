% plot the Force terms ep based on the interface
close all

if MRDGF_Eplot == 1
    MRDG_F_Ext = MRDG_F_ExtA;
elseif MRDGF_Eplot == 2
    MRDG_F_Ext = MRDG_F_ExtB;
elseif MRDGF_Eplot == 3
    MRDG_F_Ext = MRDG_F_ExtC;
elseif MRDGF_Eplot == 4
    MRDG_F_Ext = MRDG_F_ExtD;
elseif MRDGF_Eplot == 5
    MRDG_F_Ext = MRDG_F_ExtE;
elseif MRDGF_Eplot == 6
    MRDG_F_Ext = MRDG_F_ExtF;
end
%% Interface elements on the Left
facesL = length(find(MatTypeTable(2,:) == 10));
numDGfacL = facesL+nummatCG;
MRDG_F_EL = zeros(6,facesL*bCrys);

for j = nummatCG+1:numDGfacL %1:ma % 
    elemL = find(RegionOnElement (:,1) == j);
    for i = 1:length(elemL)
        MRDG_F_EL(1:6,elemL(i)-numelCG) = MRDG_F_Ext(1:6,(j));
    end
end
elemL = elemL(i);
%% Interface elements on the right 
numDGfacR = ma+facesL;
MRDG_F_ER = zeros(6,facesL*bCrys);

for j = nummat+1:numDGfacR %1:ma %
    elemR = find(RegionOnElement (:,1) == j-(nummat-nummatCG));
    for i = 1:length(elemR)
        MRDG_F_ER(1:6,elemR(i)-numelCG) = MRDG_F_Ext(1:6,(j));
    end
end

%% Interface elements on the boundary
facesB = length(find(MatTypeTable(2,:) == 12));
numDGfacB = numDGfacL+facesB;
MRDG_F_EB = zeros(6,facesB*bCrys);

for j = numDGfacL+1:numDGfacB %1:ma %
    elemB = find(RegionOnElement (:,1) == j);
    for i = 1:length(elemB)
        MRDG_F_EB(1:6,length(MRDG_F_ER)+elemB(i)-elemL) = MRDG_F_Ext(1:6,(j));
    end
end

 Fint_1L = MRDG_F_EL(1,:);
 Fint_2L = MRDG_F_EL(2,:);
 Fint_3L = MRDG_F_EL(3,:);  
 Fint_4L = MRDG_F_EL(4,:);
 Fint_5L = MRDG_F_EL(5,:);
 Fint_6L = MRDG_F_EL(6,:);
 
 Fint_1R = MRDG_F_ER(1,:);
 Fint_2R = MRDG_F_ER(2,:);
 Fint_3R = MRDG_F_ER(3,:);  
 Fint_4R = MRDG_F_ER(4,:);
 Fint_5R = MRDG_F_ER(5,:);
 Fint_6R = MRDG_F_ER(6,:);
 
 Fint_1B = MRDG_F_EB(1,:);
 Fint_2B = MRDG_F_EB(2,:);
 Fint_3B = MRDG_F_EB(3,:);  
 Fint_4B = MRDG_F_EB(4,:);
 Fint_5B = MRDG_F_EB(5,:);
 Fint_6B = MRDG_F_EB(6,:);

%  Z = Fint_1L*0;
offset = 0.07;
          
for ii_inter = 1:facesL*bCrys
    Fint_color_1L = (Fint_1L(ii_inter));%-Fint_1min)./Fint_1maxmin;
    Fint_color_2L = (Fint_2L(ii_inter));%-Fint_2min)./Fint_2maxmin;
    Fint_color_3L = (Fint_3L(ii_inter));%-Fint_3min)./Fint_3maxmin;
    Fint_color_4L = (Fint_4L(ii_inter));%-Fint_4min)./Fint_4maxmin;
    Fint_color_5L = (Fint_5L(ii_inter));%-Fint_5min)./Fint_5maxmin;
    Fint_color_6L = (Fint_6L(ii_inter));%-Fint_6min)./Fint_6maxmin;
    
    Fint_color_1R = (Fint_1R(ii_inter));%-Fint_1min)./Fint_1maxmin;
    Fint_color_2R = (Fint_2R(ii_inter));%-Fint_2min)./Fint_2maxmin;
    Fint_color_3R = (Fint_3R(ii_inter));%-Fint_3min)./Fint_3maxmin;
    Fint_color_4R = (Fint_4R(ii_inter));%-Fint_4min)./Fint_4maxmin;
    Fint_color_5R = (Fint_5R(ii_inter));%-Fint_5min)./Fint_5maxmin;
    Fint_color_6R = (Fint_6R(ii_inter));%-Fint_6min)./Fint_6maxmin;

    %     Fint_color_norm = (Fint_norm(ii_inter));%-Fint_norm_min)./Fint_norm_maxmin;
    x1 = Coordinates(NodesOnElement(numelCG+ii_inter,1),1)';
    x2 = Coordinates(NodesOnElement(numelCG+ii_inter,2),1)';
    y1 = Coordinates(NodesOnElement(numelCG+ii_inter,1),2)';
    y2 = Coordinates(NodesOnElement(numelCG+ii_inter,2),2)';
    x = [x1;x2];
    y = [y1;y2];
    z = zeros(2,2);
    
    F_1L = [Fint_color_1L;Fint_color_1L];
    F_1R = [Fint_color_1R;Fint_color_1R];

    figure(1) % Fint_1
    hold on
    surface([x,x],[y,y],z,[F_1L,F_1L],'facecol','no','edgecol','interp','linew',2);
    surface([x+offset,x+offset],[y+offset,y+offset],z,[F_1R,F_1R],'facecol','no','edgecol','interp','linew',2);
    hold off
    
    figure(2)  %Fint_2
    hold on
    F_2L = [Fint_color_2L;Fint_color_2L];
    F_2R = [Fint_color_2R;Fint_color_2R];

    surface([x,x],[y,y],z,[F_2L,F_2L],'facecol','no','edgecol','interp','linew',2);%     hold off
    surface([x+offset,x+offset],[y+offset,y+offset],z,[F_1R,F_1R],'facecol','no','edgecol','interp','linew',2);
    hold off
    
    figure(3)   %Fint_3
    hold on
    F_3L = [Fint_color_3L;Fint_color_3L];
    F_3R = [Fint_color_3R;Fint_color_3R];
    surface([x,x],[y,y],z,[F_3L,F_3L],'facecol','no','edgecol','interp','linew',2);%     hold off
    surface([x+offset,x+offset],[y+offset,y+offset],z,[F_1R,F_1R],'facecol','no','edgecol','interp','linew',2);
    hold off
    
    figure(4)   %Fint_4
    hold on
    F_4L = [Fint_color_4L;Fint_color_4L];
    F_4R = [Fint_color_4R;Fint_color_4R];
    surface([x,x],[y,y],z,[F_4L,F_4L],'facecol','no','edgecol','interp','linew',2);%     hold off
    surface([x+offset,x+offset],[y+offset,y+offset],z,[F_1R,F_1R],'facecol','no','edgecol','interp','linew',2);
    hold off
    
    figure(5)   %Fint_4
    hold on
    F_5L = [Fint_color_5L;Fint_color_5L];
    F_5R = [Fint_color_5R;Fint_color_5R];
    surface([x,x],[y,y],z,[F_5L,F_5L],'facecol','no','edgecol','interp','linew',2);%     hold off
    surface([x+offset,x+offset],[y+offset,y+offset],z,[F_1R,F_1R],'facecol','no','edgecol','interp','linew',2);
    hold off
    
    figure(6)   %Fint_4
    hold on
    F_6L = [Fint_color_6L;Fint_color_6L];
    F_6R = [Fint_color_6R;Fint_color_6R];
    surface([x,x],[y,y],z,[F_6L,F_6L],'facecol','no','edgecol','interp','linew',2);%     hold off
    surface([x+offset,x+offset],[y+offset,y+offset],z,[F_1R,F_1R],'facecol','no','edgecol','interp','linew',2);
    hold off
%     figure(7)   %ep_norm
%     hold on
%     F_n = [Fint_color_norm;Fint_color_norm];
%     surface([x,x],[y,y],z,[F_n,F_n],'facecol','no','edgecol','interp','linew',2);%     hold off
%     hold off
end

for ii_inter = facesL*bCrys+1:(facesL+facesB)*bCrys
    Fint_color_1B = (Fint_1B(ii_inter));%-Fint_1min)./Fint_1maxmin;
    Fint_color_2B = (Fint_2B(ii_inter));%-Fint_2min)./Fint_2maxmin;
    Fint_color_3B = (Fint_3B(ii_inter));%-Fint_3min)./Fint_3maxmin;
    Fint_color_4B = (Fint_4B(ii_inter));%-Fint_4min)./Fint_4maxmin;
    Fint_color_5B = (Fint_5B(ii_inter));%-Fint_5min)./Fint_5maxmin;
    Fint_color_6B = (Fint_6B(ii_inter));%-Fint_6min)./Fint_6maxmin;


    %     Fint_color_norm = (Fint_norm(ii_inter));%-Fint_norm_min)./Fint_norm_maxmin;
    x1 = Coordinates(NodesOnElement(numelCG+ii_inter,1),1)';
    x2 = Coordinates(NodesOnElement(numelCG+ii_inter,2),1)';
    y1 = Coordinates(NodesOnElement(numelCG+ii_inter,1),2)';
    y2 = Coordinates(NodesOnElement(numelCG+ii_inter,2),2)';
    x = [x1;x2];
    y = [y1;y2];
    z = zeros(2,2);
    F_1B = [Fint_color_1B;Fint_color_1B];

    figure(1) % Fint_1
    hold on
    surface([x,x],[y,y],z,[F_1B,F_1B],'facecol','no','edgecol','interp','linew',2);
    hold off
    
    figure(2)  %Fint_2
    hold on
    F_2B = [Fint_color_2B;Fint_color_2B];
    surface([x,x],[y,y],z,[F_2B,F_2B],'facecol','no','edgecol','interp','linew',2);
    hold off
    
    figure(3)   %Fint_3
    hold on
    F_3B = [Fint_color_3B;Fint_color_3B];
    surface([x,x],[y,y],z,[F_3B,F_3B],'facecol','no','edgecol','interp','linew',2);
    hold off
    
    figure(4)   %Fint_4
    hold on
    F_4B = [Fint_color_4B;Fint_color_4B];
    surface([x,x],[y,y],z,[F_4B,F_4B],'facecol','no','edgecol','interp','linew',2);
    hold off
    
    figure(5)   %Fint_4
    hold on
    F_5B = [Fint_color_5B;Fint_color_5B];
    surface([x,x],[y,y],z,[F_5B,F_5B],'facecol','no','edgecol','interp','linew',2);
    hold off
    
    figure(6)   %Fint_4
    hold on
    F_6B = [Fint_color_6B;Fint_color_6B];
    surface([x,x],[y,y],z,[F_6B,F_6B],'facecol','no','edgecol','interp','linew',2);
    hold off
%     figure(7)   %ep_norm
%     hold on
%     F_n = [Fint_color_norm;Fint_color_norm];
%     surface([x,x],[y,y],z,[F_n,F_n],'facecol','no','edgecol','interp','linew',2);%     hold off
%     hold off
end
