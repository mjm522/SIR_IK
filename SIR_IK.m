%%%%%%%%%%%%%%%%%%%
%Author: Michael Mathew, 23-03-2016
%Originally written by Rachit Sapra, this code was modified to OOP format
%%%%%%%%%%%%%%%%%%%
classdef SIR_IK
    properties
        th1Max; th2Max; th3Max; %three link manipulator max angles
        th1; th2; th3; %array for each theta_i
        a; %a_i values of DH table
        d; %d_i values of DH table
        alf; %alpha_i values of DH table
        solutionGrid; %grid for solution 
        probGrid; %probability of each grid of solution
        th1gridLen; th2gridLen; th3gridLen; %grid length of each theta_i
        targetPos; %target position;
    end
    
    properties(Constant)
        MAXDIV = 4; %number of grids = 4x4x4 = 64
    end
    
    methods
        %constructor
        function obj = SIR_IK(manInfo) %assumed that the range of angles for each DOF i is 0-theta_i_MAX
            obj.th1Max = manInfo.theta1Max; obj.th2Max = manInfo.theta2Max; obj.th3Max = manInfo.theta3Max;
            obj.a = [manInfo.a1; manInfo.a2; manInfo.a3]; %a values of i DOF
            obj.d = [manInfo.d1; manInfo.d2; manInfo.d3]; % d values of i DOF
            obj.alf = [manInfo.alf1; manInfo.alf2; manInfo.alf3].*pi/180; %alpha values of i DOF
            obj.targetPos = [manInfo.target;1];
            obj.th1gridLen = floor(obj.th1Max/(2*SIR_IK.MAXDIV)); %length of each grid = 8 for MAXDIV =4
            obj.th2gridLen = floor(obj.th2Max/(2*SIR_IK.MAXDIV));
            obj.th3gridLen = floor(obj.th3Max/(2*SIR_IK.MAXDIV));
            for i = 1:SIR_IK.MAXDIV %4 for MAXDIV = 4
                obj.th1(i) = floor(obj.th1Max*(i/SIR_IK.MAXDIV)) - obj.th1gridLen; %midpoint of each grid
                obj.th2(i) = floor(obj.th2Max*(i/SIR_IK.MAXDIV)) - obj.th2gridLen;
                obj.th3(i) = floor(obj.th3Max*(i/SIR_IK.MAXDIV)) - obj.th3gridLen;
            end
            obj.solutionGrid = []; %the solution grid, this is the grid that is resampled from time to time
            obj.probGrid = []; %probability of each solution
            for k = 1:size(obj.th3,2)
                for j = 1:size(obj.th2,2)
                    for i = 1:size(obj.th1,2)
                          errVector = obj.findError([obj.th1(i), obj.th2(j), obj.th3(k)]);
                          errVector = errVector./1000;
                          errX = errVector(1); 
                          errY = errVector(2);
                          errZ = errVector(3);
                          %after computing the error from each grid centre,
                          %lesser the error, more the probability of being a solution
                          %(- sign for value less than 1)
                          prob = exp(-errX^2*2000) *exp(-errY^2*2000) *exp(-errZ^2*2000); %combining all probabilities
                          obj.probGrid=[obj.probGrid,prob]; 
                    end
                end
            end
        end
        
        function obj = reSampler(obj) %resampling algorithm
            index = randi(SIR_IK.MAXDIV^3); %4x4x4
            obj.probGrid = obj.probGrid./sum(obj.probGrid(:)); %normalizing
            beta = 0;
            mw = max(obj.probGrid);
            for var = 1:2*SIR_IK.MAXDIV^3 %sampling is twice the grid size
                beta = beta + 2*mw.*rand(1);
                 while beta > obj.probGrid(index)
                    beta = beta-obj.probGrid(index);
                    index = index+1;
                      while (index-SIR_IK.MAXDIV^3) >= 1
                          index = index-SIR_IK.MAXDIV^3;
                      end
                end
             obj.solutionGrid = [obj.solutionGrid,index];
            end
        end
        
        function errVector = findError(obj,th)
            th = th*pi/180;
            T_total = ones(4,4); %transformation matrix from 0 to 3
%             for i = 1:3 %3 DOF find the total transformation matrix; but
%             this loop could slow down the code
%                 T = [cos(th(i)), -cos(obj.alf(i))*sin(th(i)),  sin(obj.alf(i))*sin(th(i)), obj.a(i)*cos(th(i));...
%                      sin(th(i)),  cos(obj.alf(i))*cos(th(i)), -sin(obj.alf(i))*cos(th(i)), obj.a(i)*sin(th(i));...
%                      0,           sin(obj.alf(i)),             cos(obj.alf(i)),            obj.d(i);...
%                      0,           0,                           0,                          1];
%                  T_total = T_total*T;
%             end
            %error = T_total*[0;0;0;1] - obj.targetPos;
            %errVector = error(1:3);
            errVector(1)= 70*cos(th(1))*sin(th(2)+th(3)) + 270*cos(th(1))*sin(th(2)) + 374*cos(th(1))*cos(th(2)+th(3));
            errVector(2)= 70*sin(th(1))*sin(th(2)+th(3)) + 270*sin(th(1))*sin(th(2)) + 374*sin(th(1))*cos(th(2)+th(3)) - 340.713;
            errVector(3)= -374*sin(th(2)+th(3)) + 70*cos(th(2)+th(3)) + 290 + 270*cos(th(2)) - 107.390;
        end
        
        function gridDim = findGridDim(obj,th1,th2,th3)
            gridDim.theta1_min = th1 - obj.th1gridLen;
            gridDim.theta1_max = th1 + obj.th1gridLen;
            gridDim.theta2_min = th2 - obj.th2gridLen;
            gridDim.theta2_max = th2 + obj.th2gridLen;
            gridDim.theta3_min = th3 - obj.th3gridLen;
            gridDim.theta3_max = th3 + obj.th3gridLen;
        end
        
        function gridDim = findGrid(obj,gridIndx)
           if(gridIndx > SIR_IK.MAXDIV^3) % = 64 for MAXDIV = 4
               disp('please re-run the code');
               return;
           end
           for i = 0:SIR_IK.MAXDIV^2:SIR_IK.MAXDIV^3 % 0:16:64 for MAXDIV = 4
               index3 = i/SIR_IK.MAXDIV^2+1;
            if(gridIndx > i && gridIndx <= i+SIR_IK.MAXDIV^2)
                index1 = 0; index2 = 1; 
                for j = i:SIR_IK.MAXDIV:i+SIR_IK.MAXDIV^2-SIR_IK.MAXDIV % i:4:i+12 for MAXDIV = 4
                    if(gridIndx > j && gridIndx <= j+SIR_IK.MAXDIV)
                        gridIndx = gridIndx-i;
                        gridDim = obj.findGridDim(obj.th1(gridIndx-index1),obj.th2(index2),obj.th3(index3));
                        disp(strcat('grids',num2str(j+1),'-', num2str(j+SIR_IK.MAXDIV))); %display the grid number
                    end
                    index1 = index1 + SIR_IK.MAXDIV; index2 = index2+1;
                end         
            end
                index3 = index3 + 1;
           end                  
        end
        
        function [] = parSolve(obj,gridDim)
             parfor i=1:22000 %
                theta1(i) = randi([gridDim.theta1_min gridDim.theta1_max]);
                theta2(i) = randi([gridDim.theta2_min gridDim.theta2_max]);
                theta3(i) = randi([gridDim.theta3_min gridDim.theta3_max]);
                
                errVector = obj.findError([theta1(i), theta2(i), theta3(i)]); %could be made efficient

                errX1(i) = errVector(1);
                errY1(i) = errVector(2);
                errZ1(i) = errVector(3);

                prob2(i) = exp(-(errX1(i)^2)/1000) *exp(-(errY1(i)^2)/1000) *exp(-(errZ1(i)^2)/1000); %total probability
             end
              prob2=prob2/sum(prob2(:)); %normalizing

              [~, idx]=max(prob2(:)); %find the index of the max probability
              [~, B] = ind2sub(size(prob2),idx);
              if abs(errX1(B))<=5&& abs(errY1(B))<=5 && abs(errZ1(B))<=5
                  disp('theta_1');disp(theta1(B));
                  disp('theta_2');disp(theta2(B));
                  disp('theta_3');disp(theta3(B));
              else
                  disp('No solution found, you may re-run the code');
              end
        end
        
        function obj = main(obj)
            obj = obj.reSampler(); %resample the grids
            C = unique(obj.solutionGrid);
            for m = 1:size(C,2) %which grid to redevide again
                gridIndx = C(1,m); % take one grid at a time
                gridDim = obj.findGrid(gridIndx);
                obj.parSolve(gridDim);   %find solution by multi core utilization          
            end       
        end       
    end
end

