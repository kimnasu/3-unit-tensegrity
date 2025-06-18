function grad = new_energy_gradient_plate(x, N, C_s, C_b, kb, ks, l0_s, l0_b, fixed_nodes, F_ext)

n_nodes = size(N, 2);
free_nodes = setdiff(1:n_nodes, fixed_nodes);

yield_strain = 0.25; %string's yield_strain

N_free = reshape(x, 3, []); %뉴턴 어쩌구로 계산하기 위해 x를 90x1 matrix로 변환한걸 다시
% [3x30] 배열로 바꿔놓음
N(:, free_nodes) = N_free;


%% Constraints
floor_z = min(N(3, [3, 5, 9]));  % 예: 바닥면 높이
plate_z_max = max(N(3, [26, 32, 36]));  % 초기 상단판 높이
constraint_bottom = free_nodes;
constraint_top = free_nodes;

%% Strings
u_s = N(:,C_s(:,1)) - N(:,C_s(:,2)); %each string's vector
l_s = sqrt(sum(u_s.^2,1)); %each string's length
dl_s = l_s - l0_s'; %delta_L
dE_s = zeros(3, n_nodes); %3X36 matrix


for i = 1:size(C_s,1)
    a = C_s(i,1); b = C_s(i,2);
    
    %zero dividing error control
    if l_s(i) ~= 0 
        dir = u_s(:,i)/l_s(i); %방향 벡터
    else 
        dir = zeros(3,1);
    end

    strain = dl_s(i) / l0_s(i); %epsilon = (delta_L)/(original_length)

    if strain < yield_strain
        f = ks * dl_s(i) * dir; 
    else
        f = ks * yield_strain * l0_s(i) * dir;
    end
    dE_s(:,a) = dE_s(:,a) + f; 
    dE_s(:,b) = dE_s(:,b) - f;

end

%% Bars
u_b = N(:,C_b(:,1)) - N(:,C_b(:,2)); %each bar's vector
l_b = sqrt(sum(u_b.^2,1)); %each bar's length
dl_b = l_b - l0_b'; %delta_L
dE_b = zeros(3, n_nodes);

for i = 1:size(C_b,1)
    a = C_b(i,1); b = C_b(i,2);
    if l_b(i) ~= 0
        dir = u_b(:,i)/l_b(i);
    else
        dir = zeros(3,1);
    end
    f = kb * dl_b(i) * dir;
    dE_b(:,a) = dE_b(:,a) + f;
    dE_b(:,b) = dE_b(:,b) - f;
end

%% Constraint forces using stiffness penalty
stiffness_penalty = 1e6;

% Bottom constraint (prevent from going below floor)
for i = 1:length(constraint_bottom)
    idx = constraint_bottom(i);
    z_val = N(3, idx);
    if z_val < floor_z
        dz = floor_z - z_val;
        dE_b(3, idx) = dE_b(3, idx) + stiffness_penalty * dz;
    end
end

% Top constraint (prevent from going above top plate)
for i = 1:length(constraint_top)
    idx = constraint_top(i);
    z_val = N(3, idx);
    if z_val > plate_z_max
        dz = z_val - plate_z_max;
        dE_b(3, idx) = dE_b(3, idx) - stiffness_penalty * dz;
    end
end

%% 결과: 외력 - 내부력
grad_all = dE_s + dE_b - F_ext;


% 고정된 노드는 미분 대상에서 제외
grad = reshape(grad_all(:, free_nodes), [], 1);
end
