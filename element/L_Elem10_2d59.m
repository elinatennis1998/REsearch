function [ElemFEA,ElemFEB,ElemFEC,ElemFED,ElemFEE,ElemFEF,ElemFIA,ElemFIB,ElemFIC,hr] ...
    = L_Elem10_2d59(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,...
    ndf,ndm,nst,nel,nen,MateT,isw,elem,numel,numSI,PSPS,lintt3,lintq4,lintt6,lintq9,...
    nitvms,pencoeff,BoundXmid,nummatCG,ma,MaterTypeNum,TFS,factorT,factorS,eRVE,BoundEps,BoundSig,...
    RegionsOnInterface,GrainXmid,DmatSachs,BoundSachs,DispFine,ElemE,ElemS,IBP,bCrys)
         

% DG Data Load - converts single xl,ul arrays into a left (L) and right 
% (R) side for separate treatment during calculations. Use only for DG
% elements.
CGtoDGarrays
nelLP = nelL;
nelRP = nelR;

inter = elem - (numel - numSI);
nodeAR = -1;%SurfacesI(inter,1);
nodeBR = -1;%SurfacesI(inter,2);
nodeAL = -1;%SurfacesI(inter,3);
nodeBL = -1;%SurfacesI(inter,4);

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];
        
        nelLg = (ndm+1);
        nelRg = (ndm+1);
        nstLg = ndf*nelLg;
        nstRg = ndf*nelRg;
        
        ElemKRR = zeros(nstR,nstR);
        ElemKRL = zeros(nstR,nstL);
        ElemKLR = zeros(nstL,nstR);
        ElemKLL = zeros(nstL,nstL);
        ElemFR = zeros(nstR,1);
        ElemFL = zeros(nstL,1);
       
        ElemKRRg = zeros(nstR,ndf*(ndm+1));
        ElemKRLg = zeros(nstR,ndf*(ndm+1));
        ElemKLRg = zeros(nstL,ndf*(ndm+1));
        ElemKLLg = zeros(nstL,ndf*(ndm+1));
        ElemKRgR = zeros(ndf*(ndm+1),nstR);
        ElemKRgL = zeros(ndf*(ndm+1),nstL);
        ElemKLgR = zeros(ndf*(ndm+1),nstR);
        ElemKLgL = zeros(ndf*(ndm+1),nstL);
        ElemKRgRg = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKRgLg = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKLgRg = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKLgLg = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKRgRgB = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKRgLgB = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKLgRgB = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKLgLgB = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKRgRgC = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKRgLgC = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKLgRgC = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKLgLgC = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemFRg = zeros(2*ndf*(ndm+1),1);
        ElemFLg = zeros(2*ndf*(ndm+1),1);
        ElemFEA = zeros(2*ndf*(ndm+1),1);
        ElemFEB = zeros(2*ndf*(ndm+1),1);
        ElemFEC = zeros(2*ndf*(ndm+1),1);
        ElemFED = zeros(2*ndf*(ndm+1),1);
        ElemFEE = zeros(2*ndf*(ndm+1),1);
        ElemFEF = zeros(2*ndf*(ndm+1),1);
        ElemFIA = zeros(2*ndf*(ndm+1),1);
        ElemFIB = zeros(2*ndf*(ndm+1),1);
        ElemFIC = zeros(2*ndf*(ndm+1),1);
        
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        thick = 1;
        
        xBg_mid = BoundXmid(:,ma-nummatCG); % midpoint of grain boundary; formed by InterIntegV.m
        
        % Flags for using mesoscale or microscale; these are used only in
        % the stiffness and force vectors; 0 turns off, 1 turns on
        Lmeso = mateprop(5);
        Rmeso = mateprop(6);
        Lmicr = mateprop(7);
        Rmicr = mateprop(8);
        
        %Modification of 3/21/19 (EG)
         if all(ElemFlagL==ElemFlagR)
            DiricBC = 1; % Need to flip normal vector t3 for Dirichlet edges
            fluxjc = zeros(2,1);
             fluxjx = zeros(2,1);
             fluxjy = zeros(2,1);
             jumpic = zeros(2,1);
             jumpix = zeros(2,1);
             jumpiy = zeros(2,1);
             jumpiu = zeros(2,1);
             jumpie = zeros(3,1);
        else
            DiricBC = 0;
            interI = ma-MaterTypeNum(2)+1;
%              fluxjc = CoordinatesI(2*(ndm+1)*(interI-1)+1,:)';
%              fluxjx = CoordinatesI(2*(ndm+1)*(interI-1)+2,:)';
%              fluxjy = CoordinatesI(2*(ndm+1)*(interI-1)+3,:)';
%              %             jumpic = zeros(2,1);
%              %             jumpix = zeros(2,1);
%              %             jumpiy = zeros(2,1);
%              jumpiu = CoordinatesI(2*(ndm+1)*(interI-1)+4,:)';
%              jumpie = [CoordinatesI(2*(ndm+1)*(interI-1)+5,:)'; CoordinatesI(2*(ndm+1)*(interI-1)+6,1)];
             %             jumpic = CoordinatesI(2*(ndm+1)*(interI-1)+4,:)';
             %             jumpix = CoordinatesI(2*(ndm+1)*(interI-1)+5,:)';
         end
        
         %End of modification of 3/21/19 (EG)
        
        ElemFlagLg = ElemFlag(nelL+nelR+1:nelL+nelR+ndm+1);
        ElemFlagRg = ElemFlag(nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        xlLg = xl(1:ndm,nelL+nelR+1:nelL+nelR+ndm+1);
        xlRg = xl(1:ndm,nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        ulLg = ul(1:ndf,nelL+nelR+1:nelL+nelR+ndm+1);
        ulRg = ul(1:ndf,nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        
        bf1 = 0;
        bf2 = 0;
        
        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        ElemER = matepropR(1);
        ElemvR = matepropR(2);
        lamdaR = ElemvR*ElemER/((1+ElemvR)*(1-2*ElemvR));
        muR = ElemER/(2*(1+ElemvR));
        lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
        muL = ElemEL/(2*(1+ElemvL));
        if PSPS == 'n'
        DmatR = muR*diag([2 2 1]) + lamdaR*[1; 1; 0]*[1 1 0];
        DmatL = muL*diag([2 2 1]) + lamdaL*[1; 1; 0]*[1 1 0];
        else
        DmatL = ElemEL/(1-ElemvL^2)*[1      ElemvL  0
                                  ElemvL  1      0
                                  0      0      (1-ElemvL)/2];
        DmatR = ElemER/(1-ElemvR^2)*[1      ElemvR  0
                                  ElemvR  1      0
                                  0      0      (1-ElemvR)/2];
        end
        
        NmatL = zeros(2,nstL);
        BmatL = zeros(3,nstL);
        bnAdN1 = zeros(3,nstL);
        N1 = zeros(2,nstL);
        NmatR = zeros(2,nstR);
        BmatR = zeros(3,nstR);
        bnAdN2 = zeros(3,nstR);
        N2 = zeros(2,nstR);
        
        ulresL = reshape(ulL,ndf*nelL,1);
        ulresR = reshape(ulR,ndf*nelR,1);
        
        ulresLg = reshape(ulLg,ndf*nelLg,1);
        ulresRg = reshape(ulRg,ndf*nelRg,1);
        
        % Determin bounds of integration segment
        
         InterConn2D2 % InterConn2DT %
        %End of Modification of 4/17/2019
        
        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);
        
        % Perform integration of various matrices
        
        % THIS LOOP COULD BE EVALUATED IN CLOSED-FORM
        % For separate bubble types on T and Q
        [tauL,intbL] = TauE1_2d(xlintL,DmatL,nelL,lintt6,lintq9);
        [tauR,intbR] = TauE1_2d(xlintR,DmatR,nelR,lintt6,lintq9);
        
        lint = 3;
        ib = 0;
        der = 0;
        bf = 0;
        ebL = 0;
        ebR = 0;
        intedge = 0;
        
       TaySach_case53
         
          for ie = 1:lint
         if Lmeso == 1 || Rmeso == 1
             Rmeso;
         end
         if ndm == 2
             if TFS == 2
                 DiricBC = 0;
                 fluxjc = zeros(2,1);
                 fluxjx = zeros(2,1);
                 fluxjy = zeros(2,1);
                 jumpic = zeros(2,1);
                 jumpix = zeros(2,1);
                 jumpiy = zeros(2,1);
                 jumpiu = zeros(2,1);
                 jumpie = zeros(3,1);
                 
             elseif TFS ~= 5
                 DiricBC = 0;
                 fluxjc = fluxT;
                 fluxjx = [0;0];
                 fluxjy = [0;0];
                 jumpiu = jumpSu;
                 jumpie = jumpSe; [jumpSe(3) 0];
             end
         elseif ndm == 3
             if TFS == 2
                 DiricBC = 0;
                 fluxjc = zeros(3,1);
                 fluxjx = zeros(3,1);
                 fluxjy = zeros(3,1);
                 jumpic = zeros(3,1);
                 jumpix = zeros(3,1);
                 jumpiy = zeros(3,1);
                 jumpiu = zeros(3,1);
                 jumpie = zeros(6,1);
                 
             elseif TFS ~= 5
                 DiricBC = 0;
                 fluxjc = fluxT;
                 fluxjx = [0;0;0];
                 fluxjy = [0;0;0];
                 jumpiu = jumpSu;
                 jumpie = jumpSe; [jumpSe(3) 0];
             end
             
         end
     %End of modified part EG
          end
        
%      if Rmicr + Lmicr > 0
        if nitvms == 1
            if maL == maR % Dirichlet BC assignment
                % VMS
                edgeK = tauL*ebL^2;
                gamL = eye(2);
                gamR = zeros(2,2);
                ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
            else
                % VMS
                edgeK = tauL*ebL^2 + tauR*ebR^2;
                gamL = ebL^2*(edgeK\tauL);
                gamR = ebR^2*(edgeK\tauR);
                ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
            end
        elseif nitvms == 2
        % Nitsche
        volL = getvol(xlL,nelL);
        volR = getvol(xlR,nelR);
        h = 2/(intedge/volR + intedge/volL);
        gamL = 1/2*eye(2);
        gamR = 1/2*eye(2);
        ep = pencoeff*(eye(2)*max(muL,muR)/h + (nvect'*nvect)*max(lamdaL,lamdaR)/h); %Larson CMAME version
        else
        % RFB
        tauL = RFBEfunc(xlintL,muL,lamdaL,8);
        tauR = RFBEfunc(xlintR,muR,lamdaR,8);
        edgeK = tauL + tauR;
        gamL = (edgeK\tauL);
        gamR = (edgeK\tauR);
        ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
        end

    %      end
ep = 7e3*eye(2);
% ep = ep./bCrys;
gamL = .5*eye(2);
gamR = .5*eye(2);
        
         for ie = 1:lint
            
            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
            end
            
            r = drdr*(r-roL)+eL1;
            
            if nelL == 3 || nelL == 6
                [shlL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgt(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            end
            
            rR = m*(r-eL2) + eR1;
            
            if nelR == 3 || nelR == 6
                s = 0;
            else %if nelR == 4
                s = -1;
            end
            
            if nelR == 3 || nelR == 6
                [shlR,shld,shls,be] = shlt(rR,s,nelR,nelR,der,bf);
                [QxyR, shgs, Jdet, be, xsR] = shgt(xlR(:,1:nelR),nelR,shld,shls,nen,bf,der,be);
            elseif nelR == 4 || nelR == 9
                [shlR,shld,shls,be] = shlq(rR,s,nelR,nelR,der,bf);
                [QxyR, shgs, Jdet, be, xsR] = shgq(xlR(:,1:nelR),nelR,shld,shls,nen,bf,der,be);
            end
            
            %Evaluate tangent and normal vectors
            t1 = [xsL(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
%             if DiricBC
%                 t3 = -t3;
%             end
            [tm3, tu3] = VecNormalize(t3);
            nLx = tu3(1);
            nLy = tu3(2);
            nvect = [nLx 0 nLy
                     0 nLy nLx];
            nvec = [nLx; nLy];
                    
            c1 = Wgt*tm3*drdr*thick;
            
            for i = 1:nelL
                NmatL(1,(i-1)*ndf+1) = shlL(i);
                NmatL(2,(i-1)*ndf+2) = shlL(i);
                BmatL(Bcol1,(i-1)*ndf+1) = QxyL(i,col1);
                BmatL(Bcol2,(i-1)*ndf+2) = QxyL(i,col2);
            end
            
            for i = 1:nelR
                NmatR(1,(i-1)*ndf+1) = shlR(i);
                NmatR(2,(i-1)*ndf+2) = shlR(i);
                BmatR(Bcol1,(i-1)*ndf+1) = QxyR(i,col1);
                BmatR(Bcol2,(i-1)*ndf+2) = QxyR(i,col2);
            end
            
            % Micro quantities
            bnAdN1 = gamL*nvect*DmatL*BmatL;
            bnAdN2 = gamR*nvect*DmatR*BmatR;
        
            xint = xlL*shlL; % coordinates of current integration point on edge
            xe = xint - xBg_mid; % relative coordinate to the grain edge center
            xeL = xint - xlLg(:,1); % relative coordinate to the L grain; set within GrainIntegV
            xeR = xint - xlRg(:,1); % relative coordinate to the R grain; set within GrainIntegV
            xvectL = [xeL(1) 0 xeL(2)/2
                     0 xeL(2) xeL(1)/2];
            xvectR = [xeR(1) 0 xeR(2)/2
                     0 xeR(2) xeR(1)/2];
             
         if TFS == 5 %&& DiricBC == 0
             DispFine_resR = (1-DiricBC)*reshape(DispFine(ElemFlagR,:)',ndf*nelR,1);
             DispFine_resL = reshape(DispFine(ElemFlagL,:)',ndf*nelL,1);
%              DispFine_resR = (1-DiricBC)*reshape(DispCoarse(ElemFlagR,:)',ndf*nelR,1);
%              DispFine_resL = reshape(DispCoarse(ElemFlagL,:)',ndf*nelL,1);
             jumpi =  -(NmatR*DispFine_resR - NmatL*DispFine_resL);
             
             % Flux jumps for IBP and Strong cases
             if IBP == 1
                 fluxj = -(nvect*DmatR*BmatR*DispFine_resR - nvect*DmatL*BmatL*DispFine_resL);%IBP term
                 fluxFS = zeros(2,1);
             else
                 fluxj = zeros(2,1);
                 fluxFS = (gamR*nvect*DmatR*BmatR*DispFine_resR + gamL*nvect*DmatL*BmatL*DispFine_resL);
             end
         else
                      
             tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average stress
             jumpu = NmatR*ulresR - NmatL*ulresL;        %displacement jump
             jumpi = jumpiu + (xvectR*inv(DmatR) - xvectL*inv(DmatL))*jumpie;
             %             jumpi = jumpic + xe(1)*jumpix + xe(2)*jumpiy;
             fluxj = fluxjc + xe(1)*fluxjx + xe(2)*fluxjy;
             fluxFS = zeros(2,1);
             %
             %             tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average stress
             %             jumpu = NmatR*ulresR - NmatL*ulresL;        %displacement jump

         end
%             
            % Meso quantities
            BmatLg = [0 0 1 0 0 0
                      0 0 0 1 0 0
                      0 0 0 0 1 0];
            BmatRg = BmatLg;
            bnAdN1g = gamL*nvect*DmatL*BmatLg;
            bnAdN2g = gamR*nvect*DmatR*BmatRg;
            NmatLg = [1 0 xeL(1) 0 0.5*xeL(2) 0.5*xeL(2)
                      0 1 0 xeL(2) 0.5*xeL(1) -0.5*xeL(1)];
            NmatRg = [1 0 xeR(1) 0 0.5*xeR(2) 0.5*xeR(2)
                      0 1 0 xeR(2) 0.5*xeR(1) -0.5*xeR(1)];
            tvtrg = (bnAdN1g*ulresLg + bnAdN2g*ulresRg); %average stress
            jumpug = NmatRg*ulresRg - NmatLg*ulresLg;        %displacement jump
                  
            ElemKLgLgC = ElemKLgLgC - c1*NmatLg'*bnAdN1g;
            ElemKLgRgC = ElemKLgRgC - c1*NmatLg'*bnAdN2g;
            ElemKRgLgC = ElemKRgLgC + c1*NmatRg'*bnAdN1g;
            ElemKRgRgC = ElemKRgRgC + c1*NmatRg'*bnAdN2g;
% 
            ElemKLgLg = ElemKLgLg - c1*bnAdN1g'*NmatLg;
            ElemKLgRg = ElemKLgRg + c1*bnAdN1g'*NmatRg;
            ElemKRgLg = ElemKRgLg - c1*bnAdN2g'*NmatLg;
            ElemKRgRg = ElemKRgRg + c1*bnAdN2g'*NmatRg;

            ElemKLgLgB = ElemKLgLgB + c1*(NmatLg'*ep*NmatLg);
            ElemKLgRgB = ElemKLgRgB - c1*(NmatLg'*ep*NmatRg);
            ElemKRgLgB = ElemKRgLgB - c1*(NmatRg'*ep*NmatLg);
            ElemKRgRgB = ElemKRgRgB + c1*(NmatRg'*ep*NmatRg);
            
            %Forces
            ElemFED(1:6)  = ElemFED(1:6)  + Lmeso*c1*( + bnAdN1g'*jumpi);
            ElemFED(7:12) = ElemFED(7:12) + Rmeso*c1*( + bnAdN2g'*jumpi);
            ElemFEE(1:6)  = ElemFEE(1:6)  + Lmeso*c1*( - NmatLg'*ep*jumpi);
            ElemFEE(7:12) = ElemFEE(7:12) + Rmeso*c1*( + NmatRg'*ep*jumpi);
            ElemFEF(1:6)  = ElemFEF(1:6)  - Lmeso*c1*( + NmatLg'*gamR'*fluxj);
            ElemFEF(7:12) = ElemFEF(7:12) - Rmeso*c1*( + NmatRg'*gamL'*fluxj); 
            % Version 2
            ElemFEF(1:6)  = ElemFEF(1:6)  - Lmeso*c1*( - NmatLg'*fluxFS);
            ElemFEF(7:12) = ElemFEF(7:12) - Rmeso*c1*( + NmatRg'*fluxFS);
            % Collecting fluxes:
            ElemFEC(1+2*(ie-1):ie*2) = [fluxFS(1); fluxj(1)];
            ElemFEC(7+2*(ie-1):(ie+3)*2) = [fluxFS(2); fluxj(2)];
        end

        ElemFIA = [Lmeso*Lmeso*ElemKLgLg Lmeso*Rmeso*ElemKLgRg
                   Rmeso*Lmeso*ElemKRgLg Rmeso*Rmeso*ElemKRgRg]... %Internal forces due to disp. jump
                 *[ulresLg; ulresRg];
        ElemFIB = [Lmeso*Lmeso*ElemKLgLgB Lmeso*Rmeso*ElemKLgRgB
                   Rmeso*Lmeso*ElemKRgLgB Rmeso*Rmeso*ElemKRgRgB]... %Internal forces due to tau*disp. jump
                 *[ulresLg; ulresRg];
        ElemFIC = [Lmeso*Lmeso*ElemKLgLgC Lmeso*Rmeso*ElemKLgRgC
                   Rmeso*Lmeso*ElemKRgLgC Rmeso*Rmeso*ElemKRgRgC]... %internal forces due to flux average
                 *[ulresLg; ulresRg];