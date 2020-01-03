function [ElemS,ElemE] = L_Elem12_2d53(mateprop,xl,ElemFlag,ndf,ndm,nst,nel,nen,PSPS,iprob,lintt3,lintq4,lintt6,lintq9,MateT,isw,elem,numel,numSI,nitvms,pencoeff)

        
% DG Data Load - converts single xl,ul arrays into a left (L) and right 
% (R) side for separate treatment during calculations. Use only for DG
% elements.
CGtoDGarrays

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
        
        ElemS = zeros(7,1);
        ElemE = zeros(3+9+9,1);
        t1 = zeros(3,1);
        t2 = t1;
        t3 = t1;
        thick = 1;
        
        xlLg = xl(1:ndm,nelL+nelR+1:nelL+nelR+ndm+1);
        xlRg = xl(1:ndm,nelL+nelR+ndm+1+1:nelL+nelR+ndm+1+ndm+1);
        
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
        NmatR = zeros(2,nstR);
        
        % Determin bounds of integration segment
        InterConn2Dir % InterConn2DT % 

        % Set jacobian for integration space
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
            end
        
            xint = xlL*shlL;
            xeL = xint - xlLg(:,1);
            xvectL = [xeL(1) 0 xeL(2)/2
                     0 xeL(2) xeL(1)/2];
            
            trjump = reshape( - nvect*DmatL,6,1);
            trjumpL = reshape(inv(DmatL),9,1);
            
            ElemE(1:3+9+9) = ElemE(1:3+9+9) + c1*[xint; 0; trjumpL; trjumpL];
            ElemS(1:7) = ElemS(1:7) + c1*[1; trjump];
            
        end