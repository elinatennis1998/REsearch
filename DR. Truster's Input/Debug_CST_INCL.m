% Parts of Manual Integration 
clear
%Coordinates 
elem73 = [2 2; 3 2; 2 2];
elem74 = [3 2; 4 2; 3 2];
elem76 = [4 2; 4 3; 4 2];

elem77 = [2 3; 2 4; 2 3];
elem78 = [3 4; 2 4; 3 4];
elem79 = [4 3; 4 4; 4 3];

syms x y
DirchD = [0.2*x; 0.2*y];
B_CS = [0 0 1 0 0 0
          0 0 0 1 0 0
          0 0 0 0 1 0];
delta_L = [0.341 -0.1016
           -0.1016 0.341];
delta_R = [0.659 0.1016
           0.1016 0.659];
L_73 = [-0.0058 -0.0058 0 0.0058 0.0058 0];
R_73 = [0 0.0785 0.0598 0.0598 0.0245 0.038];

L_74 = [0 0.0058 -0.0058 0.0058 -0.0058 0];
R_74 = [-0.0598 0 -0.0245 0.0598 0.0785 0.038];

L_75 = [0.0058 -0.0058 0 0 -0.0058 0.0058];
R_75 = [0.0598 0.0785 0.038 0.0598 0 0.0245];

L_76 = [0.0058 -0.0058 0 -0.0058 0 0.0058];
R_76 = [-0.0785 0 -0.0598 0.0598 -0.038 0.0245];

L_77 = [-0.0058 0.0058 0 0.0058 0 -0.0058];
R_77 = [0.0785 0.0598 0.038 0 -0.0598 -0.0245];

L_78 = [-0.0058 0 0.0058 0.0058 0 -0.0058];
R_78 = [-0.0598 -0.0598 -0.0785 0 -0.038 -0.0245];

L_79 = [-0.0058 0 0.0058 0.0058 0 0.0058];
R_79 = [-0.0598 -0.0598 -0.0785 0 -0.038 -0.0245];

L_80 = [0.0058 0 -0.0058 0.0058 -0.0058 0];
R_80 = [0 -0.0598 -0.0245 -0.0785 -0.0598 -0.038];
%% Element 29
N1 = (elem73(2,1)*elem73(3,2)-elem73(3,1)*elem73(2,2)+(elem73(2,2)-elem73(3,2))*x+(elem73(3,1)-elem73(2,1))*y);
N2 = (elem73(3,1)*elem73(1,2)-elem73(1,1)*elem73(3,2)+(elem73(3,2)-elem73(1,2))*x+(elem73(1,1)-elem73(3,1))*y);
N3 = (elem73(1,1)*elem73(2,2)-elem73(2,1)*elem73(1,2)+(elem73(1,2)-elem73(2,2))*x+(elem73(2,1)-elem73(1,1))*y);

N_73 = [N1 0 N2 0 N3 0
        0 N1 0 N2 0 N3];
n_73 = [0 0 -1
        0 -1 0];
tau_73 = [1119.2 -84.8
         -84.8 1119.2];
D_29 = 2.*[106.667 26.667 0
        26.667 106.667 0
        0 0 40];
N_CS = 1/3.*[3 0 -1.5 0 -1.5 -1.5
             0 3 0 -3 -0.75 0.75];    
%Integrate edge 73. See notes for reference
%Neighboring element is 18
D_18 = 0.5.*D_29;
%First integral of the RHS
F1_73 = vpa(int(N_CS'*tau_73*DirchD - (n_73*D_29*B_CS)'*DirchD,x,[2 3]),10); %term 1
F1_73 = vpa(subs(F1_73,y,2),10);
%Second integral of the RHS
F2_73 = vpa(int((N_CS'*tau_73*(N_73*L_73' + N_73*R_73') + ...
        (delta_L*n_73*D_29*B_CS)'*(N_73*L_73'-N_73*R_73')),x,[2 3]),10);
F2_73 = vpa(subs(F2_73,y,2),10);
F3_73 = vpa(int(N_CS'*delta_R*n_73*D_29*B_CS*L_73'-N_CS'*delta_L*n_73*D_18*B_CS*R_73',x,[2 3]),10);
F3_73 = vpa(subs(F3_73,y,2),10);
%Third integral of the RHS
F_73 = vpa(F1_73+F2_73-F3_73,10);

%Integrate edge 75. See notes for reference
%Neighboring element is 28
% tau_75 = [279.7876 -21.196
%          -21.196 279.7876];
tau_75 = [1119.2 -84.8
         -84.8 1119.2];
D_28 = 0.5.*D_29;
n_75 = [-1 0 0
        0 0 -1];  

%First integral of the RHS
F1_75 = vpa(int(N_CS'*tau_75*DirchD - (n_75*D_29*B_CS)'*DirchD,y,[2 3]),10); %term 1
F1_75 = vpa(subs(F1_75,x,2),10);
%Second integral of the RHS
F2_75 = vpa(int((N_CS'*tau_75*(N_73*L_75' + N_73*R_75') + ...
        (delta_L*n_75*D_29*B_CS)'*(N_73*L_75'-N_73*R_75')),y,[2 3]),10);  
F2_75 = vpa(subs(F2_75,x,2),10);
F3_75 = vpa(int(N_CS'*delta_R*n_75*D_29*B_CS*L_75'-N_CS'*delta_L*n_75*D_28*B_CS*R_75',y,[2 3]),10);
F3_75 = vpa(subs(F3_75,x,2),10);
%Third integral of the RHS
F_75 = vpa(F1_75+F2_75-F3_75,10);
Force_29 = F_73+F_75
%% Element 31
N1 = (elem74(2,1)*elem74(3,2)-elem74(3,1)*elem74(2,2)+(elem74(2,2)-elem74(3,2))*x+(elem74(3,1)-elem74(2,1))*y);
N2 = (elem74(3,1)*elem74(1,2)-elem74(1,1)*elem74(3,2)+(elem74(3,2)-elem74(1,2))*x+(elem74(1,1)-elem74(3,1))*y);
N3 = (elem74(1,1)*elem74(2,2)-elem74(2,1)*elem74(1,2)+(elem74(1,2)-elem74(2,2))*x+(elem74(2,1)-elem74(1,1))*y);

N_74 = [N1 0 N2 0 N3 0
        0 N1 0 N2 0 N3];
n_74 = [0 0 -1
        0 -1 0];
% tau_74 = [279.7876 21.196
%          21.196 279.7876];
tau_74 = [1119.2 84.8
         84.8 1119.2];
D_31 = 2.*[106.667 26.667 0
        26.667 106.667 0
        0 0 40];
%Integrate edge 74. See notes for reference
%Neighboring element is 20
D_20 = 0.5.*D_31;
N_20 = [ 2 - y,     0, x + y - 5,         0, 4 - x,     0
    0, 2 - y,         0, x + y - 5,     0, 4 - x];
B_20 = [0 0 1 0 -1 0
        0 -1 0 1 0 0
        -1 0 1 1 0 -1];
%First integral of the RHS
F1_74 = vpa(int(N_CS'*tau_74*DirchD - (n_74*D_31*B_CS)'*DirchD,x,[3 4]),10); %term 1
F1_74 = vpa(subs(F1_74,y,2),10);
%Second integral of the RHS
F2_74 = vpa(int((N_CS'*tau_74*(N_74*L_74' + N_74*R_74') + ...
        (delta_L*n_74*D_31*B_CS)'*(N_74*L_74'-N_74*R_74')),x,[3 4]),10);
F2_74 = vpa(subs(F2_74,y,2),10);
F3_74 = vpa(int(N_CS'*delta_R*n_74*D_31*B_CS*L_74'-N_CS'*delta_L*n_74*D_20*B_CS*R_74',x,[3 4]),10);
F3_74 = vpa(subs(F3_74,y,2),10);
%Third integral of the RHS
F_74 = vpa(F1_74+F2_74-F3_74,10);
Force_31 = F_74
%% Element 32
N1 = (elem76(2,1)*elem76(3,2)-elem76(3,1)*elem76(2,2)+(elem76(2,2)-elem76(3,2))*x+(elem76(3,1)-elem76(2,1))*y);
N2 = (elem76(3,1)*elem76(1,2)-elem76(1,1)*elem76(3,2)+(elem76(3,2)-elem76(1,2))*x+(elem76(1,1)-elem76(3,1))*y);
N3 = (elem76(1,1)*elem76(2,2)-elem76(2,1)*elem76(1,2)+(elem76(1,2)-elem76(2,2))*x+(elem76(2,1)-elem76(1,1))*y);

N_76 = [N1 0 N2 0 N3 0
        0 N1 0 N2 0 N3];
n_76 = [1 0 0
        0 0 1];
% tau_76 = [279.7876 21.196
%          21.196 279.7876];
tau_76 = [1119.2 84.8
         84.8 1119.2];
D_32 = [106.667 26.667 0
        26.667 106.667 0
        0 0 40];
    
%Integrate edge 76. See notes for reference
%Neighboring element is 33
D_33 = 2.*D_29;
N_33 = [ 7 - y - x,         0, x - 4,     0, y - 2,     0
    0, 7 - y - x,     0, x - 4,     0, y - 2];

%First integral of the RHS
F1_76 = vpa(int(N_CS'*tau_76*DirchD - (n_76*D_32*B_CS)'*DirchD,y,[2 3]),10); %term 1
F1_76 = vpa(subs(F1_76,x,4),10);
%Second integral of the RHS
F2_76= vpa(int((N_CS'*tau_76*(N_76*L_76' + N_76*R_76') + ...
        (delta_L*n_76*D_32*B_CS)'*(N_76*L_76'-N_76*R_76')),y,[2 3]),10);
F2_76 = vpa(subs(F2_76,x,4),10);
F3_76 = vpa(int(N_CS'*delta_R*n_76*D_32*B_CS*L_76'-N_CS'*delta_L*n_76*D_33*B_CS*R_76',y,[2 3]),10);
F3_76 = vpa(subs(F3_76,x,4),10);
%Third integral of the RHS
F_76 = vpa(F1_76+F2_76-F3_76,10);
Force_32 = F_76
%% Element 44
N1 = (elem79(2,1)*elem79(3,2)-elem79(3,1)*elem79(2,2)+(elem79(2,2)-elem79(3,2))*x+(elem79(3,1)-elem79(2,1))*y);
N2 = (elem79(3,1)*elem79(1,2)-elem79(1,1)*elem79(3,2)+(elem79(3,2)-elem79(1,2))*x+(elem79(1,1)-elem79(3,1))*y);
N3 = (elem79(1,1)*elem79(2,2)-elem79(2,1)*elem79(1,2)+(elem79(1,2)-elem79(2,2))*x+(elem79(2,1)-elem79(1,1))*y);

N_80 = [N1 0 N2 0 N3 0
        0 N1 0 N2 0 N3];
n_79 = [1 0 0
        0 0 1];
tau_79 = [1119.2 -84.8
         -84.8 1119.2];
D_44 = 2.*[106.667 26.667 0
        26.667 106.667 0
        0 0 40];

%Integrate edge 73. See notes for reference
%Neighboring element is 45
D_45 = 0.5.*D_44;
N_45 = [ 8 - y - x,         0, x - 4,     0, y - 3,     0
    0, 8 - y - x,     0, x - 4,     0, y - 3];

%First integral of the RHS
F1_79 = vpa(int(N_CS'*tau_79*DirchD - (n_79*D_44*B_CS)'*DirchD,y,[3 4]),10); %term 1
F1_79 = vpa(subs(F1_79,x,4),10);
%Second integral of the RHS
F2_79 = vpa(int((N_CS'*tau_79*(N_80*L_79' + N_80*R_79') + ...
        (delta_L*n_79*D_44*B_CS)'*(N_80*L_79'-N_80*R_79')),y,[3 4]),10);
F2_79 = vpa(subs(F2_79,x,4),10);
F3_79 = vpa(int(N_CS'*delta_R*n_79*D_44*B_CS*L_79'-N_CS'*delta_L*n_79*D_45*B_CS*R_79',y,[3 4]),10);
F3_79 = vpa(subs(F3_79,x,4),10);
%Third integral of the RHS
F_79 = vpa(F1_79+F2_79-F3_79,10);

%Integrate edge 80. See notes for reference
%Neighboring element is 55
% tau_80 = [279.7876 -21.196
%          -21.196 279.7876];
tau_80 = [1119.2 -84.8
         -84.8 1119.2];
D_55 = 0.5.*D_29;
n_80 = [0 0 1
        0 1 0];      
N_55 = [ 8 - y - x,         0, x - 3,     0, y - 4,     0
    0, 8 - y - x,     0, x - 3,     0, y - 4];

%First integral of the RHS
F1_80 = vpa(int(N_CS'*tau_80*DirchD - (n_80*D_44*B_CS)'*DirchD,x,[3 4]),10); %term 1
F1_80 = vpa(subs(F1_80,y,4),10);
%Second integral of the RHS
F2_80 = vpa(int((N_CS'*tau_80*(N_80*L_80' + N_80*R_80') + ...
        (delta_L*n_80*D_44*B_CS)'*(N_80*L_80'-N_80*R_80')),x,[3 4]),10);
F2_80 = vpa(subs(F2_80,y,4),10);
F3_80 = vpa(int(N_CS'*delta_R*n_80*D_44*B_CS*L_80'-N_CS'*delta_L*n_80*D_55*B_CS*R_80',x,[3 4]),10);
F3_80 = vpa(subs(F3_80,y,4),10);
%Third integral of the RHS
F_80 = vpa(F1_80+F2_80-F3_80,10);
Force_44 = F_80+F_79

%% Element 42
N1 = (elem78(2,1)*elem78(3,2)-elem78(3,1)*elem78(2,2)+(elem78(2,2)-elem78(3,2))*x+(elem78(3,1)-elem78(2,1))*y);
N2 = (elem78(3,1)*elem78(1,2)-elem78(1,1)*elem78(3,2)+(elem78(3,2)-elem78(1,2))*x+(elem78(1,1)-elem78(3,1))*y);
N3 = (elem78(1,1)*elem78(2,2)-elem78(2,1)*elem78(1,2)+(elem78(1,2)-elem78(2,2))*x+(elem78(2,1)-elem78(1,1))*y);

N_78 = [N1 0 N2 0 N3 0
        0 N1 0 N2 0 N3];
n_78 = [0 0 1
        0 1 0];
% tau_78 = [279.7876 21.196
%          21.196 279.7876];
tau_78 = [1119.2 84.8
         84.8 1119.2];     
D_42 = 2.*[106.667 26.667 0
        26.667 106.667 0
        0 0 40];

%Integrate edge 78. See notes for reference
%Neighboring element is 53
D_53 = 0.5.*D_29;

%First integral of the RHS
F1_78 = vpa(int(N_CS'*tau_78*DirchD - (n_78*D_42*B_CS)'*DirchD,x,[2 3]),10); %term 1
F1_78 = vpa(subs(F1_78,y,4),10);
%Second integral of the RHS
F2_78 = vpa(int((N_CS'*tau_78*(N_78*L_78' + N_78*R_78') + ...
        (delta_L*n_78*D_42*B_CS)'*(N_78*L_78'-N_78*R_78')),x,[2 3]),10);
F2_78 = vpa(subs(F2_78,y,4),10);
F3_78 = vpa(int(N_CS'*delta_R*n_78*D_42*B_CS*L_78'-N_CS'*delta_L*n_78*D_53*B_CS*R_78',x,[2 3]),10);
F3_78 = vpa(subs(F3_78,y,4),10);
%Third integral of the RHS
F_78 = vpa(F1_78+F2_78-F3_78,10);
Force_42 = F_78
%% Element 41
N1 = (elem77(2,1)*elem77(3,2)-elem77(3,1)*elem77(2,2)+(elem77(2,2)-elem77(3,2))*x+(elem77(3,1)-elem77(2,1))*y);
N2 = (elem77(3,1)*elem77(1,2)-elem77(1,1)*elem77(3,2)+(elem77(3,2)-elem77(1,2))*x+(elem77(1,1)-elem77(3,1))*y);
N3 = (elem77(1,1)*elem77(2,2)-elem77(2,1)*elem77(1,2)+(elem77(1,2)-elem77(2,2))*x+(elem77(2,1)-elem77(1,1))*y);

N_77 = [N1 0 N2 0 N3 0
        0 N1 0 N2 0 N3];
n_77 = [-1 0 0
        0 0 -1];
tau_77 = [1119.2 84.8
         84.8 1119.2];     
D_41 = 2.*[106.667 26.667 0
        26.667 106.667 0
        0 0 40];

%Integrate edge 77. See notes for reference
%Neighboring element is 40
D_40 = 0.5.*D_29;

%First integral of the RHS
F1_77 = vpa(int(N_CS'*tau_77*DirchD - (n_77*D_41*B_CS)'*DirchD,y,[2 3]),10); %term 1
F1_77 = vpa(subs(F1_77,x,2),10);
%Second integral of the RHS
F2_77 = vpa(int((N_CS'*tau_77*(N_77*L_77' + N_77*R_77') + ...
        (delta_L*n_77*D_41*B_CS)'*(N_77*L_77'-N_77*R_77')),y,[2 3]),10);
F2_77 = vpa(subs(F2_77,x,2),10);
F3_77 = vpa(int(N_CS'*delta_R*n_77*D_41*B_CS*L_77'-N_CS'*delta_L*n_77*D_40*B_CS*R_77',y,[2 3]),10);
F3_77 = vpa(subs(F3_77,x,2),10);
%Third integral of the RHS
F_77 = vpa(F1_77+F2_77-F3_77,10);
Force_41 = F_77
Force = vpa(Force_29+Force_31+Force_32+Force_44+Force_42+Force_41,10)   