function [ElemS,ElemE] = L_Elem12_2d51(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,PSPS,iprob,lintt3,lintq4,lintt6,lintq9)

        
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
        
        ElemS = zeros(20+6,1);
        ElemE = zeros(20,1);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        thick = 1;
        
        bf1 = 0;
        bf2 = 0;
        
        ElemFlagLg = ElemFlag(nelL+nelR+1:nelL+nelR+ndm+1);
        ElemFlagRg = ElemFlag(nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        xlLg = xl(1:ndm,nelL+nelR+1:nelL+nelR+ndm+1);
        xlRg = xl(1:ndm,nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        ulLg = ul(1:ndf,nelL+nelR+1:nelL+nelR+ndm+1);
        ulRg = ul(1:ndf,nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        
        xBg_mid = BoundXmid(:,ma-nummatCG); % midpoint of grain boundary
        
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
        
        % Determin bounds of integration segment
        InterConn2D2 % InterConn2DT % 

        % Set jacobian for integration space
        drdr = (eL2 - eL1)/drL;
        
        m = (eR2-eR1)/(eL1-eL2);
        
        % Perform integration of various matrices
        
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
                    
            if nelR == 3 || nelR == 6
                [WgtR,rR,sR] = intpntt(ie,lint,1);
                ebeR = edgebubble(rR,sR,nelR);
            elseif nelR == 4 || nelR == 9
                [WgtR,rR,sR] = intpntq(ie,lint,1);
                ebeR = edgebubbleQ(rR,sR,nelR);
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
            ebR = ebR + c1*ebeR;
            intedge = intedge + c1;
            
        end
        
        if nitvms == 1
        % VMS
        edgeK = tauL*ebL^2 + tauR*ebR^2;
        gamL = ebL^2*(edgeK\tauL);
        gamR = ebR^2*(edgeK\tauR);
        ep = pencoeff*intedge*inv(edgeK); %#ok<MINV>
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
            
            bnAdN1 = gamL*nvect*DmatL*BmatL;
            bnAdN2 = gamR*nvect*DmatR*BmatR;
        
            xint = xlL*shlL;
            xe = xint - xBg_mid;
            xeL = xint - xlLg(:,1);
            xeR = xint - xlRg(:,1);
            xvectL = [xeL(1) 0 xeL(2)/2
                     0 xeL(2) xeL(1)/2];
            xvectR = [xeR(1) 0 xeR(2)/2
                     0 xeR(2) xeR(1)/2];
            xvect = [xint(1) 0 xint(2)/2
                     0 xint(2) xint(1)/2];
            tracLR = nvect*sigLR./[volLR; volLR];
            epsV = epsLR./[volLR; volLR; volLR];
            jumpug = (disLR(1:2,2) + xvectR*epsV(1:3,2)) - (disLR(1:2,1) + xvectL*epsV(1:3,1));        %displacement jump
            
            tvtr = (bnAdN1*ulresL + bnAdN2*ulresR); %average stress
            jumpu = NmatR*ulresR - NmatL*ulresL;        %displacement jump
            
            straiLR = nvect*epsV;
            jumpsig = nvect*DmatR*BmatR*ulresR - nvect*DmatL*BmatL*ulresL;        %traction jump
            tracAve = gamL*tracLR(:,1) + gamR*tracLR(:,2);
            tracJump = tracLR(:,2) - tracLR(:,1);
            straiJump = straiLR(:,2) - straiLR(:,1);
            straiAve = gamL*straiLR(:,1) + gamR*straiLR(:,2);
            sL = DmatL*BmatL*ulresL; sR = DmatR*BmatR*ulresR;
            uL = NmatL*ulresL; uR = NmatR*ulresR;
%             HL1 = ((sL(1)-GrainSig(1,5)/9)*nLx+(sL(3)-GrainSig(3,5)/9)*nLy)*(uL(1) - GrainEps(1,5)/9*xint(1) - GrainEps(3,5)/9*xint(2));
%             HL2 = ((sL(3)-GrainSig(3,5)/9)*nLx+(sL(2)-GrainSig(2,5)/9)*nLy)*(uL(2) - GrainEps(3,5)/9*xint(1) - GrainEps(2,5)/9*xint(2));
%             HR1 = ((sR(1)-GrainSig(1,5)/9)*-nLx+(sR(3)-GrainSig(3,5)/9)*-nLy)*(uR(1) - GrainEps(1,5)/9*xint(1) - GrainEps(3,5)/9*xint(2));
%             HR2 = ((sR(3)-GrainSig(3,5)/9)*-nLx+(sR(2)-GrainSig(2,5)/9)*-nLy)*(uR(2) - GrainEps(3,5)/9*xint(1) - GrainEps(2,5)/9*xint(2));
            HL1 = ((sL(1)-GrainSig(1,maL)/9)*nLx+(sL(3)-GrainSig(3,maL)/9)*nLy)*(uL(1) - GrainEps(1,maL)/9*xint(1) - GrainEps(3,maL)/9*xint(2));
            HL2 = ((sL(3)-GrainSig(3,maL)/9)*nLx+(sL(2)-GrainSig(2,maL)/9)*nLy)*(uL(2) - GrainEps(3,maL)/9*xint(1) - GrainEps(2,maL)/9*xint(2));
            HR1 = ((sR(1)-GrainSig(1,maR)/9)*-nLx+(sR(3)-GrainSig(3,maR)/9)*-nLy)*(uR(1) - GrainEps(1,maR)/9*xint(1) - GrainEps(3,maR)/9*xint(2));
            HR2 = ((sR(3)-GrainSig(3,maR)/9)*-nLx+(sR(2)-GrainSig(2,maR)/9)*-nLy)*(uR(2) - GrainEps(3,maR)/9*xint(1) - GrainEps(2,maR)/9*xint(2));
            HillL = (nvect*DmatL*BmatL*ulresL - tracLR(:,1))'*(NmatL*ulresL - xvect*epsV(:,1));
            HillR = (-nvect*DmatR*BmatR*ulresR + tracLR(:,2))'*(NmatR*ulresR - xvect*epsV(:,2));
            HillL1 = (nvect*DmatL*BmatL*ulresL - tracLR(:,1));
            HillR1 = (-nvect*DmatR*BmatR*ulresR + tracLR(:,2));
            HillL2 = (NmatL*ulresL - xvect*epsV(:,1));
            HillR2 = (NmatR*ulresR - xvect*epsV(:,2));
            
            ElemS(1:18+8) = ElemS(1:18+8) + c1*[tvtr; (jumpu-jumpug); jumpsig; tracAve; tracJump; straiAve; straiJump; diag(gamL); diag(gamR);...
                xe(1)*(jumpu-jumpug); xe(2)*(jumpu-jumpug); xe(1)*(jumpsig-tracJump); xe(2)*(jumpsig-tracJump);];
%             ElemE(1:11) = ElemE(1:11) + c1*[1; HillL; HillR; HillL1; HillR1; HillL2; HillR2];
            ElemE(1:11) = ElemE(1:11) + c1*[1; HL1+HL2; HR1+HR2; HillL1; HillR1; HillL2; HillR2];
            
            
        end