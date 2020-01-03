% plot the Force terms for solid
close all

numfac = nummat-nummatCG;
MRDG_F_I = zeros(6,numgrain);

for j = numfac+1:nummat %1:ma % 
    elem = find(RegionOnElement (:,1) == j);
    for i = 1:length(elem)
        MRDG_F_I(1:6,elem(i)-numfac) = MRDG_F_IntA(1:6,(j));
    end
end

         Fint_1 = MRDG_F_I(1,:);
         Fint_2 = MRDG_F_I(2,:);
         Fint_3 = MRDG_F_I(3,:);  
         Fint_4 = MRDG_F_I(4,:);
         Fint_5 = MRDG_F_I(5,:);
         Fint_6 = MRDG_F_I(6,:);
         
         Fint_1max= max(Fint_1);
         Fint_1min = min(Fint_1);
         Fint_1maxmin = max(Fint_1)-min(Fint_1);
         Fint_2max= max(Fint_2);
         Fint_2min = min(Fint_2); 
         Fint_2maxmin = max(Fint_2)-min(Fint_2);  
         Fint_3max= max(Fint_3);
         Fint_3min = min(Fint_3); 
         Fint_3maxmin = max(Fint_3)-min(Fint_3);  
         Fint_4max= max(Fint_4);
         Fint_4min = min(Fint_4);  
         Fint_4maxmin = max(Fint_4)-min(Fint_4);  
         Fint_5max= max(Fint_5);
         Fint_5min = min(Fint_5);  
         Fint_5maxmin = max(Fint_5)-min(Fint_5);  
         Fint_6max= max(Fint_6);
         Fint_6min = min(Fint_6);
         Fint_6maxmin = max(Fint_6)-min(Fint_6);
         
         % Force norm
%          for ii_inter = 1:numSI %1:numel-num_locked_g   %1:numel 
%              Fint_norm(ii_inter) = sqrt(Fint_1(ii_inter)^2+Fint_2(ii_inter)^2+Fint_3(ii_inter)^2+Fint_4(ii_inter)^2+Fint_5(ii_inter)^2+Fint_6(ii_inter)^2);      
%          end

%          Fint_norm_max = max(Fint_norm);
%          Fint_norm_min = min(Fint_norm);
%          Fint_norm_maxmin = Fint_norm_max-Fint_norm_min;
         Z = Fint_1*0;
for ii_inter = 1:numSI
    Fint_color_1 = (Fint_1(ii_inter));%-Fint_1min)./Fint_1maxmin;
    Fint_color_2 = (Fint_2(ii_inter));%-Fint_2min)./Fint_2maxmin;
    Fint_color_3 = (Fint_3(ii_inter));%-Fint_3min)./Fint_3maxmin;
    Fint_color_4 = (Fint_4(ii_inter));%-Fint_4min)./Fint_4maxmin;
    Fint_color_5 = (Fint_5(ii_inter));%-Fint_5min)./Fint_5maxmin;
    Fint_color_6 = (Fint_6(ii_inter));%-Fint_6min)./Fint_6maxmin;
%     Fint_color_norm = (Fint_norm(ii_inter));%-Fint_norm_min)./Fint_norm_maxmin;
    x1 = Coordinates(NodesOnElement(numelCG+ii_inter,1),1)';
    x2 = Coordinates(NodesOnElement(numelCG+ii_inter,2),1)';
    y1 = Coordinates(NodesOnElement(numelCG+ii_inter,1),2)';
    y2 = Coordinates(NodesOnElement(numelCG+ii_inter,2),2)';
    x = [x1;x2];
    y = [y1;y2];
    z = zeros(2,2);
    F_1 = [Fint_color_1;Fint_color_1];
    figure(1) % Fint_1
    hold on
   surface([x,x],[y,y],z,[F_1,F_1],'facecol','no','edgecol','interp','linew',2);
    hold off
    figure(2)  %Fint_2
    hold on
    F_2 = [Fint_color_2;Fint_color_2];
    surface([x,x],[y,y],z,[F_2,F_2],'facecol','no','edgecol','interp','linew',2);%     hold off
    figure(3)   %Fint_3
    hold on
    F_3 = [Fint_color_3;Fint_color_3];
    surface([x,x],[y,y],z,[F_3,F_3],'facecol','no','edgecol','interp','linew',2);%     hold off
    hold off
    figure(4)   %Fint_4
    hold on
    F_4 = [Fint_color_4;Fint_color_4];
    surface([x,x],[y,y],z,[F_4,F_4],'facecol','no','edgecol','interp','linew',2);%     hold off
    hold off
    figure(5)   %Fint_4
    hold on
    F_5 = [Fint_color_5;Fint_color_5];
    surface([x,x],[y,y],z,[F_5,F_5],'facecol','no','edgecol','interp','linew',2);%     hold off
    hold off
    figure(6)   %Fint_4
    hold on
    F_6 = [Fint_color_6;Fint_color_6];
    surface([x,x],[y,y],z,[F_6,F_6],'facecol','no','edgecol','interp','linew',2);%     hold off
    hold off
%     figure(7)   %ep_norm
%     hold on
%     F_n = [Fint_color_norm;Fint_color_norm];
%     surface([x,x],[y,y],z,[F_n,F_n],'facecol','no','edgecol','interp','linew',2);%     hold off
%     hold off
end