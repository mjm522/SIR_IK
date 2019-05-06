clc; clear;
%three link manipulator DH parameters.
manInfo.theta1Max = 165; %maximum angle of theta1
manInfo.theta2Max = 110; %maximum angle of theta2
manInfo.theta3Max = 70;  %maximum angle of theta3
manInfo.alf1 = 70; %alpha 1
manInfo.alf2 = 270; %alpha 2
manInfo.alf3 = 340; %alpha 3
manInfo.d1 = 0;
manInfo.d2 = 0;
manInfo.d3 = 0;
manInfo.a1 = 0; %a1
manInfo.a2 = 0; %a2
manInfo.a3 = 0; %a3
manInfo.target = [70;70;70]; %target position [x;y;z]
myRobot = SIR_IK(manInfo);
tic;
myRobot = myRobot.main();
elap = toc %to see the time to compute solution