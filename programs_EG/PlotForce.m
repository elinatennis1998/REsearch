% % plot the stress and gap relation for the tension-shear case
% % Pinlei Chen 
% % 6/4/2014
% 
% % x axis norm of the stress plot node 6 and node 22 tension-shear case
% % node 522 and 23 for fiber -matrix case
% % node 166 and 2736 for uniform fiber-matrix
Stresstype=1;
nodeL = 3;
nodeR = 19;




normStrL = zeros(step+1,1);
Deltax = zeros(step+1,1);
Deltay = zeros(step+1,1);

stressnodeL = zeros(step+1,1);
stressnodeR = zeros(step+1,1);
step = 19;
for i = 1:step % tensioncase
   StressL = StreList(:,nodeL,i);
   stressnodeL(i+1)=StressL(Stresstype);
   normStrL(i+1) = sqrt(   StressL(1)*StressL(1) +    StressL(2)*StressL(2) ...
             +    StressL(3)*StressL(3) + 2*StressL(4)*StressL(4) ...
             + 2*StressL(5)*StressL(5) + 2*StressL(6)*StressL(6));
         
   StressR = StreList(:,nodeR,i);
   stressnodeR(i+1)=StressR(Stresstype);
   normStrR(i) = sqrt(   StressR(1)*StressR(1) +    StressR(2)*StressR(2) ...
             +    StressR(3)*StressR(3) + 2*StressR(4)*StressR(4) ...
             + 2*StressR(5)*StressR(5) + 2*StressR(6)*StressR(6));
         
   DispLx(i) = DispFine(1,nodeL,i)';
   DispRx(i) = DispFine(1,nodeR,i)';
   DispLy(i) = DispFine(2,nodeL,i)';   
   DispRy(i) = DispFine(2,nodeR,i)';   
   Deltax(i+1) = DispRx(i) - DispLx(i);
   Deltay(i+1) = DispRy(i) - DispLy(i);

%    DispLx(i) = DispList(1,nodeL,i);   
%    DispRx(i) = DispList(1,nodeR,i);
%    DispLy(i) = DispList(2,nodeL,i);   
%    DispRy(i) = DispList(2,nodeR,i);   
%    Deltax(i+1) = DispRx(i) - DispLx(i);
%    Deltay(i+1) = DispRy(i) - DispLy(i);
end
% exact solution
% DeltaS = zeros(step+1,1);
% for i = 1:step+1
%     if i ==1 ||i ==2 ||i ==3
%        DeltaS(i) = 0;
%     end
% end
ForceType=1;
nodeWForce=21;
AppliedDisp=1.5;
 
AppliedDispT=zeros(step+1,1);
AppliedDispT=linspace(0,AppliedDisp,step+1);
computedstresso1=zeros(step+1,1);
computedstresso1(2:step+1)=-squeeze(ForcList(ForceType,nodeWForce,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AppliedDispT=zeros(step+1,1);
% %  AppliedDispT=[linspace(0,0.3,13),linspace(0.375,0.33,4),linspace(0.33,1.5,47)];
%  AppliedDispT=[linspace(0,0.275,12),linspace(0.28125,0.325,8),linspace(0.35,1.5,47)];
% computedstresso1=zeros(step+1,1);
% computedstresso1(2:step+1)=-squeeze(ForcList(ForceType,nodeWForce,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 


n=51;
figure(n)
hold on
% plot(Deltax,normStrL,'ro-','LineWidth',2)
% plot(Deltax,normStrL,'b--','LineWidth',2)
% plot(Deltax,stressnodeL ,'ro-','LineWidth',2)

%      plot(AppliedDispT,computedstresso1,'r*--','LineWidth',0.5) %brittle
% %     plot(AppliedDispT,computedstresso1,'-bo','LineWidth',0.5)

%             plot(AppliedDispT,computedstresso1,'-bo','LineWidth',1,'Markersize',6)
%               plot(AppliedDispT,computedstresso1,'r*--','LineWidth',1,'Markersize',6) %brittle 
          plot(AppliedDispT,computedstresso1,'k--s','LineWidth',1,'Markersize',6)


%     %ductile

%    plot(Deltax,stressnodeL,'-bo','LineWidth',2,'Markersize',0.5) 

% plot(Deltay,normStrL,'ro-','LineWidth',2)
%plot(normStrL,Deltay,'b*-','LineWidth',2)

% hold on
% hold off
ylabel('force (N)','FontSize',14,'Fontname', 'Times New Roman')
xlabel('applied displacement (mm)','FontSize',14,'Fontname', 'Times New Roman')
box on
% grid on
% lege1=texlabel('Computed T_n vs. zeta_n');
% lege2=texlabel('Exact solution');
% h=legend(lege1,lege2,'Location','NorthEast');
% set(h,'FontSize',14,'Fontname', 'Times New Roman')
% axis([0 2.15 0 1.5])

 axis([0 1.55 0 2.5]) %example1
% axis([0 1.55 0 2.0]) %example2
set(gca,'FontSize',13,'Fontname', 'Times New Roman');
% k = 1;
% for i = 1:size(NodeTable,1)
%         if NodeTable(i,1) == 0
%           NodeBC_nodeLS(k) = i;
%           k = k+1;
%           
%         end
% end
% 
% NodeBC_nodeLS = NodeBC_nodeLS';
% for i = 1:step % x-y case 
%    Force_x = ForcList(1,NodeBC_nodeLS,i);
%    Force_y = ForcList(2,NodeBC_nodeLS,i);
%    Force_x_S(i) = sum(Force_x);
%    Force_y_S(i) = sum(Force_y);  
%    DispLx(i) = DispList(1,6,i);   
%    DispRx(i) = DispList(1,22,i);
%    DispLy(i) = DispList(2,6,i);   
%    DispRy(i) = DispList(2,22,i);   
%    Deltax(i) = DispRx(i) - DispLx(i);
%    Deltay(i) = DispRy(i) - DispLy(i);
% end
% figure()
% hold on
% plot(Deltax,Force_x_S,'ro-','LineWidth',2)
% %plot(Deltay,Force_y_S,'ro-','LineWidth',2)
% hold off
% xlabel('\delta_x (mm)','FontSize',14)
% ylabel('reaction force at x direction (N)','FontSize',14)
% box on
% grid on
% legend('F_x vs.\delta_x ','FontSize',14)
% %axis([0 0.2 0 60])
% set(gca,'FontSize',14);
% figure()
% hold on
% plot(-Deltay,-Force_y_S,'ro-','LineWidth',2)
% %plot(normStrL,Deltay,'b*-','LineWidth',2)
% hold off
% xlabel('\delta_y (mm)','FontSize',14)
% ylabel('reaction force at y direction (N)','FontSize',14)
% box on
% grid on
% legend('F_y vs.\delta_y ','FontSize',14)
% %axis([0 0.2 0 60])
% set(gca,'FontSize',14);



% % 
% 
% 
% % for i = 1:step  %fiber-matrix
% %    StressL = StreList(:,522,i);
% %    normStrL(i) = sqrt(   StressL(1)*StressL(1) +    StressL(2)*StressL(2) ...
% %              +    StressL(3)*StressL(3) + 2*StressL(4)*StressL(4) ...
% %              + 2*StressL(5)*StressL(5) + 2*StressL(6)*StressL(6));
% %    StressR = StreList(:,23,i);
% %    normStrR(i) = sqrt(   StressR(1)*StressR(1) +    StressR(2)*StressR(2) ...
% %              +    StressR(3)*StressR(3) + 2*StressR(4)*StressR(4) ...
% %              + 2*StressR(5)*StressR(5) + 2*StressR(6)*StressR(6));   
% %    DispLx(i) = DispList(1,522,i);   
% %    DispRx(i) = DispList(1,23,i);
% %    DispLy(i) = DispList(2,522,i);   
% %    DispRy(i) = DispList(2,23,i);   
% %    Deltax(i) = DispRx(i) - DispLx(i);
% %    Deltay(i) = DispRy(i) - DispLy(i);
% % end
% 
% for i = 1:step  % uniform fiber-matrix
%    StressL = StreList(:,2736,i);
%    normStrL(i) = sqrt(   StressL(1)*StressL(1) +    StressL(2)*StressL(2) ...
%              +    StressL(3)*StressL(3) + 2*StressL(4)*StressL(4) ...
%              + 2*StressL(5)*StressL(5) + 2*StressL(6)*StressL(6));
%    StressR = StreList(:,166,i);
%    normStrR(i) = sqrt(   StressR(1)*StressR(1) +    StressR(2)*StressR(2) ...
%              +    StressR(3)*StressR(3) + 2*StressR(4)*StressR(4) ...
%              + 2*StressR(5)*StressR(5) + 2*StressR(6)*StressR(6));   
%    DispLx(i) = DispList(1,2736,i);   
%    DispRx(i) = DispList(1,166,i);
%    DispLy(i) = DispList(2,2736,i);   
%    DispRy(i) = DispList(2,166,i);   
%    DispLz(i) = DispList(3,2736,i);   
%    DispRz(i) = DispList(3,166,i);   
%    Deltax(i) = DispRx(i) - DispLx(i);
%    Deltay(i) = DispRy(i) - DispLy(i);
%    Deltaz(i) = DispRz(i) - DispLz(i);   
% end
% % 
% figure()
% hold on
% plot(normStrL,DispLx,'r*-','LineWidth',2)
% plot(normStrL,DispRx,'b*-','LineWidth',2)
% plot(normStrL,DispLy,'y*-','LineWidth',2)
% plot(normStrL,DispRy,'m*-','LineWidth',2)
% hold off
% xlabel('norm of the stress \sigma')
% ylabel('Displacement')
% box on
% grid on
% legend('u_Ax','u_Bx','u_Ay','u_By')
% 
% figure()
% hold on
% plot(normStrL,Deltax,'r*-','LineWidth',2)
% plot(normStrL,Deltay,'b*-','LineWidth',2)
% hold off
% xlabel('norm of the stress \sigma')
% ylabel('\delta')
% box on
% grid on
% legend('\delta_x','\delta_y')
% 
% % figure()
% % hold on
% % plot(normStrL,DispLz,'r*-','LineWidth',2)
% % plot(normStrL,DispRz,'b*-','LineWidth',2)
% % hold off
% % xlabel('norm of the stress \sigma')
% % ylabel('Displacement')
% % box on
% % grid on
% % legend('u_Az','u_Bz')
% % 
% % figure()
% % hold on
% % plot(normStrL,Deltaz,'r*-','LineWidth',2)
% % 
% % hold off
% % xlabel('norm of the stress \sigma')
% % ylabel('\delta')
% % box on
% % grid on
% % legend('\delta_z')
% 
% for i = 1:3036
%  if NodeTable(i,1) == 1.5 &&NodeTable(i,2) == 1 && NodeTable(i,3) == 0
%   iter
%  end
% end
