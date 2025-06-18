clear; close all;
% Specify material properties
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('PLA','Poly-Flex');
material{1}='linear_elastic'; % index for material
%properties:'linear_elastic' multi elastic plastic
% cross section design coefficient
substep=100; %substep
lumped=0; % use lumped matrix 1-yes,0-no
saveimg=0; % save image or not (1) yes (0)no
savedata=0; % save data or not (1) yes (0)no
savevideo=1; % make video(1) or not(0)
gravity=0; % consider gravity 1 for yes, 0 for no
savePath=fullfile(fileparts(mfilename('fullpath')),'data_temp'); %Save
%files in same folder as this code
%% Nodal Coordinates of the structure
% Manually specify node positions of the structure.
RN=[
%% FIRST stage one and stage two
-0.7754 -0.1008 -0.023;
-0.6855  0.0155 -0.0172;
-0.7952 -0.0728 -0.0887;
-0.7059  0.0437 -0.0793;
-0.786   0.0214 -0.1057;
-0.7449 -0.0329  0.035;
-0.7414 -0.0203 -0.1339;
-0.6956 -0.0724 -0.0116;
-0.8039  0.0125 -0.0177;
-0.7142 -0.1022 -0.093;
-0.7647  0.0409 -0.0086;
-0.6704 -0.0441 -0.0849;
%% SECOND
-0.6632 -0.1026  0.0178;
-0.5765  0.0167  0.0153;
-0.6819 -0.0817 -0.0511;
-0.5972  0.0399 -0.0514;
-0.6733  0.0076 -0.0673;
-0.6345 -0.0338  0.0718;
-0.6260 -0.0281 -0.1014;
-0.5851 -0.0706  0.0301;
-0.6985 -0.0234  0.0199;
-0.5964 -0.1014 -0.0563;
-0.6569  0.0316  0.0232;
-0.5531 -0.0389 -0.0573;
%% THIRD
-0.5450 -0.1024  0.0584;
-0.4614  0.0169  0.0307;
-0.5667 -0.0865 -0.0118;
-0.4835  0.0339 -0.0330;
-0.5614  0.0104 -0.0423;
-0.5192 -0.0283  0.0957;
-0.5126 -0.0340 -0.0783;
-0.4670 -0.0681  0.0598;
-0.5852 -0.0206  0.0497;
-0.4837 -0.1013 -0.0179;
-0.5413  0.0401  0.0401;
-0.4447 -0.0505 -0.0270
]';

N=RN/2;


C_b_in = [
1 2;
3 4;
5 6;
7 8;
9 10;
11 12;
13 14;
15 16;
17 18;
19 20;
21 22;
23 24;
25 26;
27 28;
29 30;
31 32;
33 34;
35 36
];

C_s_in = [
21 13;
22 21;
18 23;
32 26;
35 26;
32 34;
33 36;
32 36;
1 6;
1 8;
1 9;
1 10;
2 6;
2 11;
2 17;
2 21;
3 5;
3 7;
3 9;
3 10;
4 5;
4 7;
4 11;
4 12;
5 9;
5 11;
6 9;
6 11;
7 10;
7 12;
8 10;
8 15;
8 21;
12 15;
12 17;
13 18;
13 20;
13 14;
18 14;
29 14;
23 14;
33 15;
19 15;
22 16;
17 16;
19 16;
23 16;
24 17;
23 18;
19 22;
19 24;
20 22;
20 27;
20 33;
24 27;
24 29;
25 30;
25 32;
25 33;
25 34;
26 27;
31 27;
34 28;
29 28;
31 28;
35 28;
36 29;
35 30;
26 30;
33 30;
35 31;
34 31;
];
% Convert the above matrices into full connectivity matrices.
C_b = tenseg_ind2C(C_b_in,N);%%
C_s = tenseg_ind2C(C_s_in,N);
C=[C_b;C_s];
[ne,nn]=size(C); % ne:No.of element;nn:No.of node

% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);

%% Boundary constraints
pinned_X=([3; 5; 9])'; pinned_Y=([3; 5; 9])'; pinned_Z=([3; 5; 9])';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);
%generate group index
gr={(1:6);(7:12);(13:18);(19:42);(43:66);(67:89)};
Gp=tenseg_str_gp(gr,C); % generate group matrix
%% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);
%SVD of equilibrium matrix
[U1,U2,V1,V2,S]=tenseg_svd1(A_1ag);
%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;
%prestress design
index_gp=1; % number of groups with designed force
t=load('t.txt');% force vector
t = t(1:90); 
t_gp=pinv(Gp)*t;
%% cross sectional design
A_b=0.00282^2*ones(18,1);
r_b=A_b/2;
A_s=((pi/4)*0.001^2)*ones(72,1);
r_s=A_s/2;

% cross sectional of all members
% Index of bar and string elements in the full connectivity matrix C
index_b = 1:size(C_b,1);                          % bar index: 1~18
index_s = size(C_b,1)+1 : size(C,1);              % string index: 19~90
I3=eye(ne);
Ind_b=I3(:,index_b); % index matrix for bar
Ind_s=I3(:,index_s); % index matrix for string
A=[Ind_b,Ind_s]*[A_b;A_s]; % cross sectional area
A_gp=pinv(Gp)*A;
radius=[Ind_b,Ind_s].*[r_b;r_s]; %radius
r_gp=pinv(Gp)*radius;
E=[Ind_b,Ind_s]*[Eb*ones(numel(index_b),1);Es*ones(numel(index_s),1)];
%Young's modulus vector
% members' force & rest length
l0=E.*A.*l./(t+E.*A);

bar_length = zeros(size(C_b,1), 1);
for i = 1:size(C_b,1)
    nodes = find(C_b(i,:) ~= 0);
    p1 = N(:, nodes(1));
    p2 = N(:, nodes(2));
    bar_length(i) = norm(p2 - p1);  % meter
end

bar_length_mm = bar_length * 1000;  % mm 단위 변환

fprintf("📏 Bar Index\tLength (mm)\n");
fprintf("-----------------------------\n");
for i = 1:length(bar_length_mm)
    fprintf("Bar %2d\t\t%.2f mm\n", i, bar_length_mm(i));
end
%t = load('t.txt');

%fprintf("\n📌 Bar 위치 정보:\n");
%fprintf("-----------------------------\n");
%for i = 1:size(C_b_in,1)
    %n1 = C_b_in(i,1);
    %n2 = C_b_in(i,2);
    %fprintf("Bar %2d: Node %2d ↔ Node %2d\n", i, n1, n2);
%end

%figure;
%hold on;

% bar 그리기 + bar 번호 표시
for i = 1:size(C_b,1)
    idx = find(C_b(i,:) ~= 0);   % 연결된 노드 인덱스
    p1 = N(:, idx(1));
    p2 = N(:, idx(2));
    mid = (p1 + p2) / 2;         % 중심점

    % bar 그리기 (검정색)
    plot3([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], 'k', 'LineWidth', 2);

    % bar 번호 표시 (파란 글씨)
    text(mid(1), mid(2), mid(3), sprintf('Bar %d', i), ...
        'FontSize', 10, 'Color', 'blue', 'FontWeight', 'bold');
end

axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
view(3);





