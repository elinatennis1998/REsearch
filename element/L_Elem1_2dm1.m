function ElemF = L_Elem1_3dm1(mateprop,nodeA,nodeB,elem,edge,traction,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen)

        
        
        ElemF = zeros(nst,1);
        
        % Reorder nodes for corresponding face of integration
        SurfOrientEtoS

        % Determine bounds of integration
        
        if nel == 4 || nel == 9
            
        dr = 2;
        ro = -1;
        
        % Upper Limit
%         if nodeA == ElemFlag(2)
            eR2 = 1;
%         else %
%             eR2 = 0;
%         end
        % Lower Limit
%         if nodeB == ElemFlag(1)
            eR1 = -1;
%         else %nodeA == ElemFlag(5)
%             eR1 = 0;
%         end
        
        elseif nel == 3 || nel == 6
            
        dr = 1;
        ro = 0;
            
        % Upper Limit
%         if nodeA == ElemFlag(2)
            eR2 = 1;
%         else %nodeA == ElemFlagR(5)
%             eR2 = 1/2;
%         end
        % Lower Limit
%         if nodeB == ElemFlag(1)
            eR1 = 0;
%         else %nodeA == ElemFlag(5)
%             eR1 = 1/2;
%         end
        
        end
        
        % Set jacobian for integration space
        drdr = (eR2 - eR1)/dr;
        
        lint = 4;
%         % Load Gauss Points for quadrature
%         if enrich == 1
%             [rlist, rWgts, rnum] = GaussPoints(pr+2);
%             slist = -1;
%             sWgts = 1;
%             snum = 1;
%         else
%             [rlist, rWgts, rnum] = GaussPoints(pr+1);
%             slist = -1;
%             sWgts = 1;
%             snum = 1;
%         end

%         lamda = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
%         mu = PatchE/(2*(1+Patchv));
%         thick = 1;
        ideriv = 0;
        der = 0;
        bf = 0;

        Nmat = zeros(2,ndf*nel);
        
        for je = 1:lint

            if nel == 3 || nel == 6
                [Wgt,r,s] = intpntt(je,lint,1);
            else %if nel == 4
                [Wgt,r,s] = intpntq(je,lint,1);
            end
            r = drdr*(r-ro)+eR1;

            if nel == 3 || nel == 6
                [shl,shld,shls,be] = shlt(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgt(xl,nel,shld,shls,nen,bf,der,be);
            elseif nel == 4 || nel == 9
                [shl,shld,shls,be] = shlq(r,s,nel,nel,der,bf);
                [shg, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
            end
                    
            %Evaluate tangent and normal vectors
            t1 = [xs(:,1); 0];
            [tm1, tu1] = VecNormalize(t1);
            t2 = [0; 0; 1];
            tm2 = 1;
            tu2 = t2';
            t3 = VecCrossProd(t1,t2);
            [tm3, tu3] = VecNormalize(t3);
            Traction = zeros(1,2);
            
            if edge < 0 % negative edge denotes pressure
                pressur = -traction(1); % pressure is defined positive in compression
                Traction = pressur*tu3(1:2);
            else
                Traction = traction(1:2);
            end
            
            c1 = Wgt*tm3*drdr*thick;
            
            for ie = 1:nel
                
              Nmat(1,(ie-1)*ndf+1) = shl(ie);
              Nmat(2,(ie-1)*ndf+2) = shl(ie);
                 
            end
            
            ElemF = ElemF + c1*Nmat'*Traction';

        end %ie
        
        % Reorder nodes back to the orientation of the element
        SurfOrientStoE
        ElemF(1:ndf*nel) = ElemF(ilist2);