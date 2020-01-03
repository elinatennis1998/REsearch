function [ElemK,ElemF,hr] = L_Elem11_2d100(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,PSPS,iprob)
         
PatchE = mateprop(1);
Patchv = mateprop(2);
thick = mateprop(3);
areaG = mateprop(4);

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];



        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);
%         Nmat = zeros(2,ndf*nel);
%         Bmat = zeros(3,ndf*nel);

        % Load Gauss Points for quadrature
        lint = 1;
 
        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        
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
        Bmat = [0 0 1 0 0 0
                0 0 0 1 0 0
                0 0 0 0 1 0];
        c1 = areaG*thick;
        
            
            sigma = Dmat*(Bmat*ulres(1:nel*ndf));
            
            ElemF = ElemF - c1*Bmat'*sigma;
            ElemK = ElemK + c1*Bmat'*Dmat*Bmat;
            
ElemK;