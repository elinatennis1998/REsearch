% Two dimenisional plane stress/strain element, coarse scale


%% Global declarations
% Place here any commands that should be executed
% whenever the element routine is called.
if ~exist('PSPS','var')
PSPS = 'n'; %plane stresS/straiN flag
end
if ~exist('iprob','var')
iprob = 0;
end


%% Task Switch - Perform various element functions
switch isw 
%%    
    case 1 % Setup up elemental degrees of freedom (OPTIONAL)
        
% Purpose: Rearrange order of dofs on element and/or flag certain
%          dofs to be presecribed away
% Called by: pmatin.m
        
        if ndf > 2
            
            for i = 3:ndf
                lie(i,1) = 0;
            end
            
        end
        
        istv = 7; % number of stresses per node for post-processing
        iste = 7; % number of stresses per element

%%
    case {3,6,21} % Stiffness and internal force vector (REQUIRED)
        
        if exist('IBP','var') && IBP == 1
            [ElemK,ElemF,hr] = L_Elem11_2d03BConly(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,PSPS,iprob);
        else
            [ElemK,ElemF,hr] = L_Elem11_2d03(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,PSPS,iprob);
        end

%%        
    case 26 % Element Stress (OPTIONAL)
        
        ElemS = L_Elem11_2d26(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,nestr,PSPS,iprob);
        
%%
    case 59 % MRDG forces
        
        [ElemFE,ElemFI,hr] = L_Elem11_2d59(mateprop,ul,xl,ElemFlag,hr,nh1,nh2,nh3,ndf,ndm,nst,nel,nen,PSPS,iprob);
%         MRDG_F_ExtA(1:6,ma) = ElemFE;
        MRDG_F_IntA(1:6,ma) = ElemFI;
        
        
end %Task Switch


PatchE = mateprop(1);
Patchv = mateprop(2);
thick = mateprop(3);
areaG = mateprop(4);

Bcol1 = [1; 3];
Bcol2 = [2; 3];
col1 = [1; 2];
col2 = [2; 1];

switch isw %Task Switch
    
    case 1
%%
    case 3
%%
    case 26
        
        ElemS = zeros(nestr,1);
        stresID = [1 2 3 0 1 3];

        lam = Patchv*PatchE/((1+Patchv)*(1-2*Patchv));
        mu = PatchE/(2*(1+Patchv));
        thick = 1;
        fbx = 0;
        fby = 0;
        fbz = 0;
        if PSPS == 'n'
        Dmat = mu*diag([2 2 1]) + lam*[1; 1; 0]*[1 1 0];
        else
        Dmat = PatchE/(1-Patchv^2)*[1      Patchv  0
                                  Patchv  1      0
                                  0      0      (1-Patchv)/2];
        end
        Bmat = [0 0 1 0 0 0
                0 0 0 1 0 0
                0 0 0 0 1 0];
        I1 = [1; 1; 0];
        
        % Load Guass Integration Points

            lint = 1;
            nint = 1;
        
        der = 0;
        bf = 0;
        ib = 0;
        
        ulres = reshape(ul,nst,1);

        %Stress Loop
        ll = 1;
            
            epsil = Bmat*ulres(1:ndf*nel);
            sigma = Dmat*epsil;
            
            for stres = 1:npstr
            
            if stres <= 3 % stress components
                sigmas = sigma(stresID(stres));
            elseif stres >= 5
                if stres <= 6 % principal stresses
                    if PSPS == 'n'
                        sigz = lam*(epsil(1)+epsil(2));
                    else
                        sigz = 0;
                    end
                    sigma2 = [sigma(1) sigma(3) 0; sigma(3) sigma(2) 0; 0 0 sigz];
                    psig = eig(sigma2);
                    sigmas = psig(stresID(stres));
                else % hydrostatic stress
                    if PSPS == 'n'
                        sigz = lam*(epsil(1)+epsil(2));
                    else
                        sigz = 0;
                    end
                    sigmas = 1/3*(sigma'*I1 + sigz);
                end
            else % von Mises stress
                if PSPS == 'n'
                    sigz = lam*(epsil(1)+epsil(2));
                else
                    sigz = 0;
                end
                trs = sigma'*I1 + sigz;
                dsig = sigma - 1/3*trs*I1;
                sigmas = sqrt(3/2*(dsig'*dsig));
            end
            
            ElemS(stres) = sigmas;
            
            end

        
end %Task Switch
