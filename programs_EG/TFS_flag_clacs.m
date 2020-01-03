%Elina Geut
%Script for jump calculations for TFS flags 
%Created 8/18/2019
%Part of case 53 and TaylorSet & SachsSet

%Part from TaylorSet (EG Modification)
         if TFS == 1
             if ndm == 2
                 trjump = reshape(trjump(1:6),2,3);
             else
                 trjump = reshape(trjump(1:36),6,6);
             end

             fluxT = factorT*trjump*eRVE';
            
             %Part from TaylorSet end
             
             %Part from Sachs set
             if ndm == 2
                 jumpSu = zeros(2,1);
                 jumpSe = zeros(3,1);
             elseif ndm == 3
                 jumpSu = zeros(3,1);
                 jumpSe = zeros(6,1);
             end
         elseif TFS == 3
             BoundEps(1:length(ElemE),interI) = BoundEps(1:length(ElemE),interI) + ElemE;
             BoundSig(1:length(ElemS),interI) = BoundSig(1:length(ElemS),interI) + ElemS;
             if ndm == 2
                 BoundSachs = BoundEps(4:21,interI)./(ones(18,1)*BoundSig(1,interI));
                 matLR = RegionsOnInterface(interI,2:3);
                 trjumpi = BoundSachs(1:18,1);
                 iDmatL = reshape(trjumpi(1:9),3,3);
                 iDmatR = reshape(trjumpi(10:18),3,3);
                 xint = GrainXmid(1:2,matLR(1)) - GrainXmid(1:2,matLR(2));
                 xvect = [xint(1) 0 xint(2)/2
                     0 xint(2) xint(1)/2];
                 fluxT = [0;0];
             elseif ndm == 3
%                  BoundSachs = BoundEps(5:28,interI)./(ones(18,1)*BoundSig(1,interI));
%              matLR = RegionsOnInterface(interI,2:3);
%              trjumpi = BoundSachs(1:18,1);
%              iDmatL = reshape(trjumpi(1:9),3,3);
%              iDmatR = reshape(trjumpi(10:18),3,3);
             xint = GrainXmid(1:3,matLR(1)) - GrainXmid(1:3,matLR(2));
             xvect = [xint(1) 0 0 xint(2)/2 xint(3)/2
                 0 xint(2) 0 xint(1)/2  xint(3)/2
                 0 0 xint(3) xint(1)/2  xint(1)/2];
             fluxT = [0;0;0];
             end
             jumpSu = factorS*xvect*eRVE';
             jumpSe = factorS*DmatSachs*eRVE';
             
         elseif TFS == 4
              if ndm == 2
                 trjump = reshape(trjump(1:6),2,3);
             else
                 trjump = reshape(trjump(1:18),3,6);
             end
             fluxT = factorT*trjump*eRVE';
             BoundEps(1:length(ElemE),interI) = BoundEps(1:length(ElemE),interI) + ElemE;
             BoundSig(1:length(ElemS),interI) = BoundSig(1:length(ElemS),interI) + ElemS;
             if ndm == 2
                 BoundSachs = BoundEps(4:21,interI)./(ones(18,1)*BoundSig(1,interI));
                 matLR = RegionsOnInterface(interI,2:3);
                 trjumpi = BoundSachs(1:18,1);
                 iDmatL = reshape(trjumpi(1:9),3,3);
                 iDmatR = reshape(trjumpi(10:18),3,3);
                 xint = GrainXmid(1:2,matLR(1)) - GrainXmid(1:2,matLR(2));
                 xvect = [xint(1) 0 xint(2)/2
                     0 xint(2) xint(1)/2];
             elseif ndm == 3
%              BoundSachs = BoundEps(4:21,interI)./(ones(18,1)*BoundSig(1,interI));
%              matLR = RegionsOnInterface(interI,2:3);
%              trjumpi = BoundSachs(1:18,1);
%              iDmatL = reshape(trjumpi(1:9),3,3);
%              iDmatR = reshape(trjumpi(10:18),3,3);
             xint = GrainXmid(1:3,matLR(1)) - GrainXmid(1:3,matLR(2));
             xvect = [xint(1) 0 0 xint(2)/2 xint(3)/2
                       0 xint(2) 0 xint(1)/2 xint(3)/2
                       0 0 xint(3) xint(1)/2 xint(2)/2];
             end
             jumpSu = factorS*xvect*eRVE';
             jumpSe = factorS*DmatSachs*eRVE';         
         end
         %End part from Sachs model