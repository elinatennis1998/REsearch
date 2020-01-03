% plot the Force terms ep based on the interface
close all

numDGfac = ma-num_locked_g;
MRDG_F_I = zeros(6,numSI);

for j = nummatCG+1:numDGfac
    elem = find(RegionOnElement (:,1) == j);
    for i = 1:count
        MRDG_F_I(1,elem-numelCG) = MRDG_F_Int(1,(j));
        MRDG_F_I(2,elem-numelCG) = MRDG_F_Int(2,(j));
        MRDG_F_I(3,elem-numelCG) = MRDG_F_Int(3,(j));
        MRDG_F_I(4,elem-numelCG) = MRDG_F_Int(4,(j));
        MRDG_F_I(5,elem-numelCG) = MRDG_F_Int(5,(j));
        MRDG_F_I(6,elem-numelCG) = MRDG_F_Int(6,(j));
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
         for ii_inter = 1:numSI
             Fint_norm(ii_inter) = sqrt(Fint_1(ii_inter)^2+Fint_2(ii_inter)^2+Fint_3(ii_inter)^2+Fint_4(ii_inter)^2+Fint_5(ii_inter)^2+Fint_6(ii_inter)^2);      
         end

         Fint_norm_max = max(Fint_norm);
         Fint_norm_min = min(Fint_norm);
         Fint_norm_maxmin = Fint_norm_max-Fint_norm_min;
for ii_inter=1:numSI
    Fint_color_1 = (Fint_1(ii_inter)-Fint_1min)./Fint_1maxmin;
    Fint_color_2 = (Fint_2(ii_inter)-Fint_2min)./Fint_2maxmin;
    Fint_color_3 = (Fint_3(ii_inter)-Fint_3min)./Fint_3maxmin;
    Fint_color_4 = (Fint_4(ii_inter)-Fint_4min)./Fint_4maxmin;
    Fint_color_5 = (Fint_5(ii_inter)-Fint_5min)./Fint_5maxmin;
    Fint_color_6 = (Fint_6(ii_inter)-Fint_1min)./Fint_6maxmin;
    Fint_color_norm = (Fint_norm(ii_inter)-Fint_norm_min)./Fint_norm_maxmin;        
    figure(1) % Fint_1
    hold on
    plot([Coordinates(NodesOnElement(ii_inter,1),1),Coordinates(NodesOnElement(ii_inter,2),1)],[Coordinates(NodesOnElement(ii_inter,1),2),Coordinates(NodesOnElement(ii_inter,2),2)],'Color',[Fint_color_1,0,0],'LineWidth',4)
    hold off
    figure(2)  %Fint_2
    hold on
    plot([Coordinates(NodesOnElement(ii_inter,1),1),Coordinates(NodesOnElement(ii_inter,2),1)],[Coordinates(NodesOnElement(ii_inter,1),2),Coordinates(NodesOnElement(ii_inter,2),2)],'Color',[Fint_color_2,0,0],'LineWidth',4)
    hold off
    figure(3)   %Fint_3
    hold on
    plot([Coordinates(NodesOnElement(ii_inter,1),1),Coordinates(NodesOnElement(ii_inter,2),1)],[Coordinates(NodesOnElement(ii_inter,1),2),Coordinates(NodesOnElement(ii_inter,2),2)],'Color',[Fint_color_3,0,0],'LineWidth',4)
    hold off
    figure(4)   %Fint_4
    hold on
    plot([Coordinates(NodesOnElement(ii_inter,1),1),Coordinates(NodesOnElement(ii_inter,2),1)],[Coordinates(NodesOnElement(ii_inter,1),2),Coordinates(NodesOnElement(ii_inter,2),2)],'Color',[Fint_color_4,0,0],'LineWidth',4)
    hold off
    figure(5)   %Fint_4
    hold on
    plot([Coordinates(NodesOnElement(ii_inter,1),1),Coordinates(NodesOnElement(ii_inter,2),1)],[Coordinates(NodesOnElement(ii_inter,1),2),Coordinates(NodesOnElement(ii_inter,2),2)],'Color',[Fint_color_5,0,0],'LineWidth',4)
    hold off
    figure(6)   %Fint_4
    hold on
    plot([Coordinates(NodesOnElement(ii_inter,1),1),Coordinates(NodesOnElement(ii_inter,2),1)],[Coordinates(NodesOnElement(ii_inter,1),2),Coordinates(NodesOnElement(ii_inter,2),2)],'Color',[Fint_color_6,0,0],'LineWidth',4)
    hold off
    figure(7)   %ep_norm
    hold on
    plot([Coordinates(NodesOnElement(ii_inter,1),1),Coordinates(NodesOnElement(ii_inter,2),1)],[Coordinates(NodesOnElement(ii_inter,1),2),Coordinates(NodesOnElement(ii_inter,2),2)],'Color',[Fint_color_norm,0,0],'LineWidth',4)
    hold off
end