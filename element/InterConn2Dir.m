% Tim Truster
% 09/07/2019
%

        
        truncsectors = 1;

        if nelL == 3 || nelL == 6
        xlintL = zeros(2,3);
        nelLB = 3;
        else
        xlintL = zeros(2,4);
        nelLB = 4;
        end
        % Determine bounds of integration, left
        
        if nelLB == 4
            
            t1 = [(xlL(1,2)-xlL(1,1)); (xlL(2,2)-xlL(2,1)); 0];
            t2 = [(xlL(1,4)-xlL(1,1)); (xlL(2,4)-xlL(2,1)); 0];
            t3 = VecCrossProd(t1,t2);
            aL = VecNormalize(t3);
            hL = sqrt((xlL(1,2)-xlL(1,1))^2+(xlL(2,2)-xlL(2,1))^2)/aL;
            
            drL = 2;
            roL = -1;

            % Upper Limit
%             if nodeAL == ElemFlagL(1) || nodeAL == 0
                eL1 = -1;
                xlintL(:,1) = xlL(:,1);
                xlintL(:,4) = xlL(:,4);
%             end
            % Lower Limit
%             if nodeBL == ElemFlagL(2) || nodeBL == 0
                eL2 = 1;
                xlintL(:,2) = xlL(:,2);
                xlintL(:,3) = xlL(:,3);
%             end
            
            if truncsectors
            % shorten support for bubble
            lside23 = norm(xlintL(:,2)-xlintL(:,3),2);
            lside14 = norm(xlintL(:,1)-xlintL(:,4),2);
            lside12 = norm(xlintL(:,2)-xlintL(:,1),2);
            shortyn = (lside12 < lside23) && (lside12 < lside14);
            if shortyn && lside14 < lside23
                lintratio = lside12/lside14;
                lintxy1 = (xlintL(:,4)-xlintL(:,1))*lintratio + xlintL(:,1);
                POUxi = POU_Coord(lintxy1(1),lintxy1(2),xlintL,0,nelLB);
%                 shl = shlq(-1,POUxi(2),nelL,nelL,0,0);
%                 lintxy1 = xlintL*shl;
                xlintL(:,4) = lintxy1;
                shl = shlq(1,POUxi(2),nelLB,nelLB,0,0);
                lintxy2 = xlintL*shl;
                xlintL(:,3) = lintxy2;
            elseif shortyn && lside14 >= lside23
                lintratio = lside12/lside23;
                lintxy1 = (xlintL(:,3)-xlintL(:,2))*lintratio + xlintL(:,2);
                POUxi = POU_Coord(lintxy1(1),lintxy1(2),xlintL,0,nelLB);
%                 shl = shlq(1,POUxi(2),nelL,nelL,0,0);
%                 lintxy1 = xlintL*shl;
                xlintL(:,3) = lintxy1;
                shl = shlq(-1,POUxi(2),nelLB,nelLB,0,0);
                lintxy2 = xlintL*shl;
                xlintL(:,4) = lintxy2;
            end
            end
        
        elseif nelLB == 3
            
            t1 = [(xlL(1,2)-xlL(1,1)); (xlL(2,2)-xlL(2,1)); 0];
            t2 = [(xlL(1,3)-xlL(1,1)); (xlL(2,3)-xlL(2,1)); 0];
            t3 = VecCrossProd(t1,t2);
            aL = VecNormalize(t3);
            aL = aL/2;
            hL = sqrt((xlL(1,2)-xlL(1,1))^2+(xlL(2,2)-xlL(2,1))^2)/aL;
            
            drL = 1;
            roL = 0;

            % Upper Limit
%             if nodeAL == ElemFlagL(1) || nodeAL == 0
                eL1 = 0;
                xlintL(:,1) = xlL(:,1);
%             end
            % Lower Limit
%             if nodeBL == ElemFlagL(2) || nodeBL == 0
                eL2 = 1;
                xlintL(:,2) = xlL(:,2);
%             end
            xlintL(:,3) = xlL(:,3);
            
            if truncsectors
            % shorten support for bubble
            lside23 = norm(xlintL(:,2)-xlintL(:,3),2);
            lside13 = norm(xlintL(:,1)-xlintL(:,3),2);
            lside12 = norm(xlintL(:,2)-xlintL(:,1),2);
            shortyn = (lside12 < lside23) && (lside12 < lside13);
            if shortyn && lside13 < lside23
                lintratio = lside12/lside13;
                lintxy1 = (xlintL(:,3)-xlintL(:,1))*lintratio + xlintL(:,1);
                shl = shlt(1-lintratio,lintratio,nelLB,nelLB,0,0);
                lintxy2 = xlintL*shl;
                xlintL(:,3) = (lintxy1 + lintxy2)/2;
            elseif shortyn && lside13 >= lside23
                lintratio = lside12/lside23;
                lintxy1 = (xlintL(:,3)-xlintL(:,2))*lintratio + xlintL(:,2);
                POUxi = POU_Coord(lintxy1(1),lintxy1(2),xlintL,1,nelLB);
                shl = shlt(0,POUxi(2),nelLB,nelLB,0,0);
                lintxy2 = xlintL*shl;
                xlintL(:,3) = (lintxy1 + lintxy2)/2;
            end
            end
        
        end
        