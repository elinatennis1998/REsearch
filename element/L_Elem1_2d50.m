function [ElemS,ElemE,ElemD,ElemV,ElemW] = L_Elem1_2d50(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,PSPS,iprob,lintt3,lintq4,lintt6,lintq9)

PatchE = mateprop(1);
Patchv = mateprop(2);
thick = mateprop(3);

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

        ElemS = zeros(3,1);
        ElemE = zeros(3,1);
        ElemD = zeros(2,1);
        ElemV = zeros(2,1);
        ElemW = zeros(1,1);
        Nmat = zeros(2,ndf*nel);
        Bmat = zeros(4,ndf*nel);

        % Load Gauss Points for quadrature
        if nel == 3
            lint = lintt3;%13;
        elseif nel == 4
            lint = lintq4;
        elseif nel == 6
            lint = lintt6;%13;
        elseif nel == 9
            lint = lintq9;
        end
 
        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        
        cvec = zeros(3,1); %temperature vector, see Hughes book
        if PSPS == 'n'
        Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        else
        Dmat = PatchE/(1-Patchv^2)*[1      Patchv  0
                                  Patchv  1      0
                                  0      0      (1-Patchv)/2];
        end

        der = 0;
        bf = 0;
        ib = 0;
        
        ulres = reshape(ul,ndf*nen,1);
        
        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,litr,lits] = intpntt(je,lint,ib);
                [shl,shld,shls,be] = shlt(litr,lits,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [Wgt,litr,lits] = intpntq(je,lint,ib);
                [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
                    
            c1 = Wgt*Jdet*thick;
            
            % Form B matrix
            for ie = 1:nel
                
              Nmat(1,(ie-1)*ndf+1) = shl(ie);
              Nmat(2,(ie-1)*ndf+2) = shl(ie);
                
              Bmat(Bcol1,(ie-1)*ndf+1) = shg(ie,col1);
              Bmat(Bcol2,(ie-1)*ndf+2) = shg(ie,col2);
              Bmat(4,(ie-1)*ndf+1:(ie-1)*ndf+2) = [Bmat(3,(ie-1)*ndf+1) -Bmat(3,(ie-1)*ndf+2)];
                 
            end
            
            xint = xl(:,1:nel)*shl;
            
            displac = Nmat*ulres(1:nel*ndf);
            epsil = (Bmat*ulres(1:nel*ndf));
            sigma = Dmat*epsil(1:3);
            
            %Modification for engineering strain 9/5/2019
            epsil(3) = 0.5*epsil(3);
            rotation = 0.5.*epsil(4);
            ElemD = ElemD + c1*displac;
            ElemE = ElemE + c1*epsil(1:3);
            ElemS = ElemS + c1*[sigma];
            ElemV = ElemV + c1*[1; sigma'*epsil(1:3)];
            %Rotation Integration added 9/5/2019
            ElemW = ElemW + c1*rotation;
            
        end %je