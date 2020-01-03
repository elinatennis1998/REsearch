function [ElemK,ElemF,hr] = L_Elem12_2d100(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,...
    ndf,ndm,nst,nel,nen,MateT,isw,elem,numel,numSI,PSPS,lintt3,lintq4,lintt6,lintq9,...
    nitvms,pencoeff,BoundXmid,nummatCG,ma,MaterTypeNum,TFS,factorT,factorS,eRVE,BoundEps,BoundSig,...
    RegionsOnInterface,GrainXmid,DmatSachs,BoundSachs)
         

% DG Data Load - converts single xl,ul arrays into a left (L) and right 
% (R) side for separate treatment during calculations. Use only for DG
% elements.

maL = mateprop(1);
nelL = mateprop(3);
nstL = nelL*ndf;
ElemFlagL = ElemFlag(1:nelL);
xlL = xl(1:ndm,1:nelL);
if isw > 1 && isw ~= 53
ulL = ul(1:ndf,1:nelL);
end
matepropL = MateT(maL,:);

nelLP = nelL;



Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];
        
        nelLg = (ndm+1);
        nelRg = (ndm+1);
        nstLg = ndf*nelLg;
        nstRg = ndf*nelRg;
        
        ElemK = zeros(nst,nst);
        ElemF = zeros(nst,1);
        ElemKLL = zeros(nstL,nstL);
        ElemFL = zeros(nstL,1);
       
        ElemKLRg = zeros(nstL,ndf*(ndm+1));
        ElemKLLg = zeros(nstL,ndf*(ndm+1));
        ElemKRgL = zeros(ndf*(ndm+1),nstL);
        ElemKLgL = zeros(ndf*(ndm+1),nstL);
        ElemKRgRg = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKRgLg = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKLgRg = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemKLgLg = zeros(ndf*(ndm+1),ndf*(ndm+1));
        ElemFRg = zeros(ndf*(ndm+1),1);
        ElemFLg = zeros(ndf*(ndm+1),1);
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
        
         %End of modification of 3/21/19 (EG)
        
        ElemFlagLg = ElemFlag(nelL+1:nelL+ndm+1);
        ElemFlagRg = ElemFlag(nelL+ndm+1+1:nelL+ndm+1+ndm+1);
        xlLg = xl(1:ndm,nelL+1:nelL+ndm+1);
        xlRg = xl(1:ndm,nelL+ndm+1+1:nelL+ndm+1+ndm+1);
        ulLg = ul(1:ndf,nelL+1:nelL+ndm+1);
        ulRg = ul(1:ndf,nelL+ndm+1+1:nelL+ndm+1+ndm+1);
        
        bf1 = 0;
        bf2 = 0;
        
        ElemEL = matepropL(1);
        ElemvL = matepropL(2);
        lamdaL = ElemvL*ElemEL/((1+ElemvL)*(1-2*ElemvL));
        muL = ElemEL/(2*(1+ElemvL));
        if PSPS == 'n'
        DmatL = muL*diag([2 2 1]) + lamdaL*[1; 1; 0]*[1 1 0];
        else
        DmatL = ElemEL/(1-ElemvL^2)*[1      ElemvL  0
                                  ElemvL  1      0
                                  0      0      (1-ElemvL)/2];
                                              
        end
        
        NmatL = zeros(2,nstL);
        BmatL = zeros(3,nstL);
        bnAdN1 = zeros(3,nstL);
        N1 = zeros(2,nstL);
        
        ulresL = reshape(ulL,ndf*nelL,1);
        
        ulresLg = reshape(ulLg,ndf*nelLg,1);
        ulresRg = reshape(ulRg,ndf*nelRg,1);
        
        %% Determine bounds of integration segment
        InterConn2Dir
        
        %% Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        % Perform integration of various matrices
        
        % For separate bubble types on T and Q
        [tauL,intbL] = TauE1_2d(xlintL,DmatL,nelL,lintt6,lintq9);
        
        lint = 3;
        ib = 0;
        der = 0;
        bf = 0;
        ebL = 0;
        ebR = 0;
        intedge = 0;
        
        % THIS LOOP COULD BE EVALUATED IN CLOSED-FORM
        for ie = 1:lint
            
% For separate bubble types on T and Q
            if nelL == 3 || nelL == 6
                [Wgt,r,s] = intpntt(ie,lint,1);
                ebeL = edgebubble(r,s,nelL);
            elseif nelL == 4 || nelL == 9
                [Wgt,r,s] = intpntq(ie,lint,1);
                ebeL = edgebubbleQ(r,s,nelL);
            end
            
            r = drdr*(r-roL)+eL1;
            
            if nelL == 3 || nelL == 6
                [shlL,shld,shls,be] = shlt(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgt(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            elseif nelL == 4 || nelL == 9
                [shlL,shld,shls,be] = shlq(r,s,nelL,nelL,der,bf);
                [QxyL, shgs, Jdet, be, xsL] = shgq(xlL(:,1:nelL),nelL,shld,shls,nen,bf,der,be);
            end
                    
            %Evaluate tangent and normal vectors
            t1 = [xsL(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            nLx = tu3(1);
            nLy = tu3(2);
            nvect = [nLx nLy];
                    
            c1 = Wgt*tm3*drdr*thick;
            
            ebL = ebL + c1*ebeL;
            intedge = intedge + c1;
        end
        
        % Hack for now
            nvect = [nLx 0 nLy
                     0 nLy nLx];
            trjump = 0*reshape(- nvect*DmatL,6,1);
            trjumpL = 0*reshape(inv(DmatL),9,1);
        
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
            
            %Elina Geut
             for i = 1:nelL
                BmatL(1,(i-1)*ndf+1) = shld(i);
                BmatL(2,(i-1)*ndf+2) = shld(i);
             end
            %End
            
          end         
         
          for ie = 1:lint
         if Lmeso == 1 || Rmeso == 1
             Rmeso;
         end
         if ndm == 2
                 DiricBC = 1; % Need to flip normal vector t3 for Dirichlet edges
                 fluxjc = zeros(2,1);
                 fluxjx = zeros(2,1);
                 fluxjy = zeros(2,1);
                 jumpic = zeros(2,1);
                 jumpix = zeros(2,1);
                 jumpiy = zeros(2,1);
                 jumpiu = zeros(2,1);
                 jumpie = zeros(3,1);
         elseif ndm == 3
                 DiricBC = 1; % Need to flip normal vector t3 for Dirichlet edges
                 fluxjc = zeros(3,1);
                 fluxjx = zeros(3,1);
                 fluxjy = zeros(3,1);
                 jumpic = zeros(3,1);
                 jumpix = zeros(3,1);
                 jumpiy = zeros(3,1);
                 jumpiu = zeros(3,1);
                 jumpie = zeros(6,1);
             
         end
     %End of modified part EG
          end
        
        edgeK = tauL*ebL^2;
        gamL = eye(2);
        gamR = zeros(2,2);
        ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
        
        %% main integration loop
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
            
            % Micro quantities
            bnAdN1 = gamL*nvect*DmatL*BmatL;
        
            xint = xlL*shlL; % coordinates of current integration point on edge
            xe = xint - xBg_mid; % relative coordinate to the grain edge center
            xeL = xint - xlLg(:,1); % relative coordinate to the L grain; set within GrainIntegV
            xeR = xint - xlRg(:,1); % relative coordinate to the R grain; set within GrainIntegV
            xvectL = [xeL(1) 0 xeL(2)/2
                     0 xeL(2) xeL(1)/2];
            xvectR = [xeR(1) 0 xeR(2)/2
                     0 xeR(2) xeR(1)/2];
            load('DispFine','DispFine');         
          if TFS == 5 %&& DiricBC == 0
             DispFine_resL = reshape(DispFine(ElemFlagL,:)',ndf*nelL,1);
%              DispFine_resL = reshape(DispFine(ElemFlagL,:)',ndf*nelL,1);
             jumpi =  -(- NmatL*DispFine_resL);
             fluxj = zeros(2,1);
             fluxFS = (gamL*nvect*DmatL*BmatL*DispFine_resL);
             
         else
                      
             tvtr = (bnAdN1*ulresL); %average stress
             jumpu = - NmatL*ulresL;        %displacement jump
             jumpi = jumpiu + (- xvectL*inv(DmatL))*jumpie;
             %             jumpi = jumpic + xe(1)*jumpix + xe(2)*jumpiy;
             fluxj = fluxjc + xe(1)*fluxjx + xe(2)*fluxjy;
             fluxFS = zeros(2,1);
             %
             %             tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average stress
             %             jumpu = NmatR*ulresR - NmatL*ulresL;        %displacement jump

         end
         
            ElemKLL = ElemKLL - c1*NmatL'*bnAdN1;
% 
            ElemKLL = ElemKLL - c1*bnAdN1'*NmatL;

            ElemKLL = ElemKLL + c1*(NmatL'*ep*NmatL);
            
            % add interface prescribed flux
%             ElemFL = ElemFL - c1*( + NmatL'*gamR'*fluxj);           
            ElemFL = ElemFL + c1*( - NmatL'*ep*jumpi);
            ElemFL = ElemFL + c1*( + bnAdN1'*jumpi);
            
            % Meso quantities
            BmatLg = [0 0 1 0 0 0
                      0 0 0 1 0 0
                      0 0 0 0 1 0];
            BmatRg = BmatLg;
            bnAdN1g = gamL*nvect*DmatL*BmatLg;
            NmatLg = [1 0 xeL(1) 0 0.5*xeL(2) 0.5*xeL(2)
                      0 1 0 xeL(2) 0.5*xeL(1) -0.5*xeL(1)];
            NmatRg = [1 0 xeR(1) 0 0.5*xeR(2) 0.5*xeR(2)
                      0 1 0 xeR(2) 0.5*xeR(1) -0.5*xeR(1)];
%             tvtrg = (bnAdN1g*ulresLg + bnAdN2g*ulresRg); %average stress
%             jumpug = NmatRg*ulresRg - NmatLg*ulresLg;        %displacement jump
                  
            ElemKLgLg = ElemKLgLg - c1*NmatLg'*bnAdN1g;
            ElemKRgLg = ElemKRgLg + c1*NmatRg'*bnAdN1g;
% 
            ElemKLgLg = ElemKLgLg - c1*bnAdN1g'*NmatLg;
            ElemKLgRg = ElemKLgRg + c1*bnAdN1g'*NmatRg;

            ElemKLgLg = ElemKLgLg + c1*(NmatLg'*ep*NmatLg);
            ElemKLgRg = ElemKLgRg - c1*(NmatLg'*ep*NmatRg);
            ElemKRgLg = ElemKRgLg - c1*(NmatRg'*ep*NmatLg);
            ElemKRgRg = ElemKRgRg + c1*(NmatRg'*ep*NmatRg);
            
            % Meso-micro quantities

            ElemKLgL = ElemKLgL - c1*NmatLg'*bnAdN1;
            ElemKRgL = ElemKRgL + c1*NmatRg'*bnAdN1;
% 
            ElemKLgL = ElemKLgL - c1*bnAdN1g'*NmatL;
%             ElemKRgL = ElemKRgL - c1*bnAdN2g'*NmatL;

            ElemKLgL = ElemKLgL + c1*(NmatLg'*ep*NmatL);
            ElemKRgL = ElemKRgL - c1*(NmatRg'*ep*NmatL);
            
            % Micro-meso quantities

            ElemKLLg = ElemKLLg - c1*NmatL'*bnAdN1g;
%             ElemKLRg = ElemKLRg - c1*NmatL'*bnAdN2g;
% 
            ElemKLLg = ElemKLLg - c1*bnAdN1'*NmatLg;
            ElemKLRg = ElemKLRg + c1*bnAdN1'*NmatRg;

            ElemKLLg = ElemKLLg + c1*(NmatL'*ep*NmatLg);
            ElemKLRg = ElemKLRg - c1*(NmatL'*ep*NmatRg);
            
            %Forces
%             ElemFLg = ElemFLg - c1*( + NmatLg'*gamR'*fluxj);
            ElemFRg = ElemFRg - c1*( + NmatRg'*gamL'*fluxj); 
            ElemFLg = ElemFLg - c1*( - NmatLg'*fluxFS);
            ElemFRg = ElemFRg - c1*( + NmatRg'*fluxFS);
            ElemFLg = ElemFLg + c1*( - NmatLg'*ep*jumpi);
            ElemFRg = ElemFRg + c1*( + NmatRg'*ep*jumpi);
            ElemFLg = ElemFLg + c1*( + bnAdN1g'*jumpi);
%             ElemFRg = ElemFRg + c1*( + bnAdN2g'*jumpi);
            
        end

            ElemK = [Lmicr*Lmicr* ElemKLL Lmicr*Lmeso* ElemKLLg Lmicr*Rmeso* ElemKLRg
                     Lmeso*Lmicr*ElemKLgL Lmeso*Lmeso*ElemKLgLg Lmeso*Rmeso*ElemKLgRg
                     Rmeso*Lmicr*ElemKRgL Rmeso*Lmeso*ElemKRgLg Rmeso*Rmeso*ElemKRgRg];
            ElemF = [Lmicr*ElemFL; Lmeso*ElemFLg; Rmeso*ElemFRg]-ElemK*[ulresL; ulresLg; ulresRg];
if elem == 81
    elem;
end
        ElemK;