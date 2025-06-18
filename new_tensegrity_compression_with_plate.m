% tensegrity_compression_with_plate.m
clear; clc;

%% Material properties
Eb = 3.6e9;
Es = 0.1e9;
Ab = 7.9524e-6;
As = 0.7854e-6;

%% Load node coordinates and connectivity
RN = [
-0.134600591871953  -0.140586680378839  -0.757662765330182;
-0.116842210651314  -0.0198074697366638 -0.675575237817924;
-0.0727466230078326 -0.117656504600978  -0.791438845247083;
-0.0586214402319784  0.00351103470616943 -0.70925608166615;
-0.0495923390045865 -0.0243229140116465 -0.791438845247083;
-0.182028725590841  -0.0678481900346446 -0.720725735747515;
-0.0153090293547824 -0.0653823596410833 -0.750670394159766;
-0.128694536493285  -0.107598669342485  -0.679040442560313;
-0.139726279727272  -0.0287379711799528 -0.791438845247083;
-0.0541149530809364 -0.143184108991829  -0.711157762233714;
-0.139531458546271   0.00208382659811394 -0.753020245240896;
-0.050556205298223  -0.0825855644867857 -0.670233298668637;
-0.152638131560937  -0.134300344312972  -0.639799076516216;
-0.127229478130217  -0.0111948176770089 -0.562630043591605;
-0.0877889932746916 -0.118587837963493  -0.672683972863157;
-0.0648123099145373  0.00684450932010562 -0.597188686409671;
-0.0657804020951007 -0.0301203657792593 -0.672819175013265;
-0.196457200625369  -0.0609976981094333 -0.605451714959826;
-0.024890140018025  -0.0654223907018532 -0.630910422900698;
-0.147767223997532  -0.097753613373069  -0.562886719140336;
-0.157665805485418  -0.0569799268260224 -0.678778757563564;
-0.0669091507172746 -0.134275628342843  -0.588759569929379;
-0.149992651698636   0.000126561677908221 -0.64074690667888;
-0.0543172124675264 -0.0698689913790545 -0.550349047853661;
-0.169203098385202  -0.1257329042166    -0.516219786809947;
-0.119724065702916  -0.00431436290138193 -0.4421;
-0.103918863858659  -0.115233680018213  -0.551989536788321;
-0.0608288596948995  0.00765671489178731 -0.481924978129671;
-0.0681742083051282 -0.0202239623193319 -0.558579160612062;
-0.196972110065699  -0.0483068506670683 -0.488258595454187;
-0.0255632682021851 -0.0642396889494656 -0.515039252564791;
-0.153562607156382  -0.087561056314449  -0.4421;
-0.164485658910204  -0.0467217891393379 -0.562241080349431;
-0.0824034586926709 -0.126216154764827  -0.470985076485303;
-0.14344386308956    0.01539745696511    -0.524793756118488;
-0.0633065557475773 -0.0741803887028924 -0.4421;
]';
N = RN / 2;

plate_nodes_target = [26, 32, 36]; %변위를 바꿀 노드

C_b= [ 
          1 2;3 4;5 6;7 8;9 10;11 12; 
          13 14;15 16;17 18;19 20;21 22;23 24; 
          25 26;27 28;29 30;31 32;33 34;35 36 
          ];

C_s= [      
1 6;1 8;1 9;1 10; 2 6;2 11;2 17;2 21;3 5;3 7;3 9;3 10;4 5;4 7;4 11;4 12; 
5 9;5 11;6 9;6 11;7 10;7 12;8 10;8 15;8 21;12 15;12 17;13 18;13 20;13 21;13 22; 
14 18;14 29;14 23;14 33;15 19;15 22;16 17;16 19;16 23;16 24;17 23;18 21;18 23; 
19 22;19 24;20 22;20 27;20 33;24 27;24 29;25 30;25 32;25 33;25 34;26 32;26 35;26 36; 
27 31;27 34;28 29;28 31;28 35;28 36;29 35;30 26;30 33;30 35;31 34;31 36;32 34;32 36 
];


%% 고정 노드 설정
fixed_nodes = [3, 5, 9, 26, 32, 36];
free_nodes = setdiff(1:size(N,2), fixed_nodes);
x0 = reshape(N(:,free_nodes), [], 1); %free node들의 위치를 1차원 벡터로 생성


%% Constarints : 모든 노드들이 바닥면 높이와 상단판의 높이를 초과하지 않도록 제한

floor_z = min(N(3, [3, 5, 9]));  % 예: 바닥면 높이
plate_z_max = max(N(3, [26, 32, 36]));  % 초기 상단판 높이

constraint_bottom = free_nodes;
constraint_top = free_nodes;

%% Rest length and stiffness
l0_b = vecnorm(N(:,C_b(:,1)) - N(:,C_b(:,2)), 2, 1)';
l0_s = 0.95 * vecnorm(N(:,C_s(:,1)) - N(:,C_s(:,2)), 2, 1)';
ks = mean((Es * As) ./ l0_s);
kb = mean((Eb * Ab) ./ l0_b);

%% fsolve로 초기 equilibrium 해 구하기
options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-5,'StepTolerance',1e-6,'MaxFunctionEvaluations',1e5);
x_sol = fsolve(@(x) new_energy_gradient_plate(x, N, C_s, C_b, kb, ks, l0_s, l0_b, fixed_nodes, zeros(3, size(N,2))), x0, options);
%x_sol : prestressed로 인해 변한 displacement matrix [x; y; z] -> 90x1 matrix

%% 비디오 저장 설정
v = VideoWriter('deformation_plate.avi'); open(v);
logfile = fopen('fsolve_plate_log.txt', 'w');
fprintf(logfile, 'Step\tDz\tExitFlag\tIterations\tFuncEvals\tFinalNorm\n');

%% Plate 압축 단계 수행
displacement_steps = linspace(0, -0.08, 100);

target_bars = 1:18;
bar_forces_log = zeros(length(displacement_steps), length(target_bars));

for step = 1:length(displacement_steps)
    dz = displacement_steps(step);
    N_step = N;
    N_step(3,plate_nodes_target) = N(3,plate_nodes_target) + dz;

    [x_sol, ~, exitflag, output] = fsolve(@(x) ...
        new_energy_gradient_plate(x, N_step, C_s, C_b, kb, ks, l0_s, l0_b, fixed_nodes, zeros(3,size(N,2))), ...
        x0, options);

    N_sol = N_step;
    N_sol(:,free_nodes) = reshape(x_sol, 3, []);
    x0 = reshape(N_sol(:,free_nodes), [], 1);

    %compression of bars
    for j = 1:length(target_bars)
        idx = target_bars(j);
        a = C_b(idx,1);
        b = C_b(idx,2);
        u = N_sol(:,a) - N_sol(:,b);
        l_now = norm(u); %변형된 길이
        l0 = l0_b(idx);
        kb_i = Eb * Ab / l0;  % 개별 강성
        bar_forces_log(step, j) = kb_i * (l_now - l0);  % 스칼라 힘 저장
    end

    grad_val = new_energy_gradient_plate(x_sol, N_step, C_s, C_b, kb, ks, l0_s, l0_b, fixed_nodes, zeros(3,size(N,2)));
    grad_norm = norm(grad_val);
    fprintf(logfile, '%d\t%.5f\t%d\t%d\t%d\t%.3e\n', step, dz, exitflag, output.iterations, output.funcCount, grad_norm);
    
    %% Video
    figure(1);
    clf;
    set(gcf, 'Units', 'pixels');
    screenSize = get(0, 'ScreenSize');
    set(gcf, 'Position', [100 100 screenSize(3)-200 screenSize(4)-200]);
    %set(gcf, 'Position', [100 100 1920 1080])
    for vview = 1:3
        subplot(1,3,vview); hold on; axis equal;
        for i = 1:size(C_b,1)
            plot3(N_sol(1,C_b(i,:)), N_sol(2,C_b(i,:)), N_sol(3,C_b(i,:)), 'k', 'LineWidth', 2);
        end
        for i = 1:size(C_s,1)
            plot3(N_sol(1,C_s(i,:)), N_sol(2,C_s(i,:)), N_sol(3,C_s(i,:)), 'r');
        end
        scatter3(N_sol(1,:), N_sol(2,:), N_sol(3,:), 10, 'b', 'filled');
        view_vec = {[0 0], [90 0], [90 90]};
        title_str = {'Front (Y-Z)', 'Side (X-Z)', 'Top (X-Y)'};
        view(view_vec{vview});
        title(title_str{vview});
        axis([-0.11 0.01 -0.1 0.02 -0.41 -0.2]);
        xlabel('X'); ylabel('Y'); zlabel('Z');
    end
    drawnow;
    frame_all = getframe(gcf);
    rotate3d on
    writeVideo(v, frame_all);
end

%% Compression force of Bar: plot
one_stage_bar_forces_log = bar_forces_log(:, 1:6);
two_stage_bar_forces_log = bar_forces_log(:, 7:12);
three_stage_bar_forces_log = bar_forces_log(:, 13:18);
target_bars_forces_log = bar_forces_log(:, [6, 7, 16]);
figure(2);
plot(-displacement_steps, abs(target_bars_forces_log), 'LineWidth', 1.5);  % 부호 반전 및 절댓값
legend('Bar 6', 'Bar 7', 'Bar 16');
xlabel('Z Displacement (m)');
ylabel('Compression Force Magnitude (N)');
title('Force vs. Displacement for Selected Bars');
grid on;

figure(3);
plot(-displacement_steps, abs(one_stage_bar_forces_log), 'LineWidth', 1.5);  % 부호 반전 및 절댓값
legend('Bar 1', 'Bar 2', 'Bar 3', 'Bar 4', 'Bar 5', 'Bar 6');
xlabel('Z Displacement (m)');
ylabel('Compression Force Magnitude (N)');
title('Force vs. Displacement for one stage Bars');
grid on;

figure(4);
plot(-displacement_steps, abs(two_stage_bar_forces_log), 'LineWidth', 1.5);  % 부호 반전 및 절댓값
legend('Bar 7', 'Bar 8', 'Bar 9', 'Bar 10', 'Bar 11', 'Bar 12');
xlabel('Z Displacement (m)');
ylabel('Compression Force Magnitude (N)');
title('Force vs. Displacement for two stage Bars');
grid on;

figure(5);
plot(-displacement_steps, abs(three_stage_bar_forces_log), 'LineWidth', 1.5);  % 부호 반전 및 절댓값
legend('Bar 13', 'Bar 14', 'Bar 15', 'Bar 16', 'Bar 17', 'Bar 18');
xlabel('Z Displacement (m)');
ylabel('Compression Force Magnitude (N)');
title('Force vs. Displacement for three stage Bars');
grid on;
writematrix(bar_forces_log, 'bar_forces_log.csv');
fclose(logfile); close(v);