function [ElemK,ElemF,hr] = L_Elem1_3d03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,ma,ieFEAP)
         
PatchE = mateprop(1);
Patchv = mateprop(2);

Bcol1 = [1; 4; 6];
Bcol2 = [2; 4; 5];
Bcol3 = [3; 5; 6];
col1 = [1; 2; 3];
col2 = [2; 1; 3];
col3 = [3; 2; 1];
iemat = ieFEAP(1:ndf,ma);

        ElemK = zeros(ndf*nel);
        ElemF = zeros(ndf*nel,1);

        %Set integration number
        if nel ~= 6
            if nel == 4
                lint = 1;4;11;5;16;
            elseif nel == 8
                lint = 8;1000; %1000 for body force problem 
            elseif nel == 10
                lint = 14;
            elseif nel == 18
                linttw = 13;
                lintlw = 3;
                lint = linttw*lintlw;
            else
                lint = 27;
            end
        else
            linttw = 7;
            lintlw = 2;
            lint = linttw*lintlw;
        end
        ib = 0;
        bf = 0;
        der = 0;

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        fbx = 0;
        fby = 0;
        fbz = 0;
        Dmat = mu*diag([2 2 2 1 1 1]) + lam*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
        Nmat = zeros(3,ndf*nel);
        Bmat = zeros(6,ndf*nel);
        rhspul = reshape(ul,ndf*nen,1);

        for l = 1:lint                    

            % Evaluate 1-D basis functions at integration points
            if nel == 4 || nel == 10
              [w,ss] =  int3d_t(l,lint,ib);
              [shl,shld,shls,be] = shltt(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be] = shgtt(xl,nel,shld,shls,nel,bf,der,be);
            elseif nel == 6 || nel == 18 % wedge
              [w,ss] =  intpntw(l,linttw,lintlw,ib);
              [shl,shld,shls,be] = shlw(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be,xs] = shgw(xl,nel,shld,shls,nel,bf,der,be);
            else
              [w,ss] =  intpntb(l,lint,ib);
              [shl,shld,shls,be] = shlb(ss,nel,nel,der,bf);
              [shg,shgs,Jdet,be,xs] = shgb(xl,nel,shld,shls,nel,bf,der,be);
            end
            
            if Jdet < 0
                elem
                Jdet %#ok<NOPTS>
            end
            c1 = Jdet*w*thick;

            Fvec = [fbx; fby; fbz];
            
            % Form B matrix
            for i = 1:nel
                
              Nmat(1,(i-1)*ndf+iemat(1)) = shl(i);
              Nmat(2,(i-1)*ndf+iemat(2)) = shl(i);
              Nmat(3,(i-1)*ndf+iemat(3)) = shl(i);
                
            Bmat(Bcol1,(i-1)*ndf+iemat(1)) = shg(i,col1);
            Bmat(Bcol2,(i-1)*ndf+iemat(2)) = shg(i,col2);
            Bmat(Bcol3,(i-1)*ndf+iemat(3)) = shg(i,col3);
                 
            end
            
            ElemF = ElemF + c1*(Nmat'*Fvec - Bmat'*Dmat*(Bmat*rhspul(1:ndf*nel)));
            
            ElemK = ElemK + c1*(Bmat'*Dmat*Bmat);

        end %je
ElemK;