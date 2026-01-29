clear; clc; close all;
%constants:
q = 1.602e-19; m0 = 9.109e-31; k = 1.38e-23; %C,kg,J/K
T = 300; %K
% Euler method helps to guess differentiation 
% with tiny steps(dt)
%new val = old val + pace*time step
dt = 1e-15; %time step
t = 0 : dt : 12e-12; %in pikosec
N = length(t); %total steps

%v(t + dt) = v(t) + dv/dt * dt;
%w(t + dt) = w(t) + dw/dt * dt;

%initially,
v = zeros(1, N); %m/s
v(1) = 0;

w = zeros(1, N); %J --- w0= 3/2kt when 300k
w_init_eV = 0.038; %38 meV bdry
w(1) = w_init_eV * q; %in J


%E- field outside:
E_field = zeros(1, N);
for i = 1:N
    t_ps = t(i) * 1e12; %piko to num
switch true
    case t_ps >= 0 && t_ps <= 2
        E_val = 8;
    case (t_ps >= 2 & t_ps <= 3) ...
            | (t_ps >= 5 & t_ps <= 10)
        E_val = 1;
    otherwise % (t_ps >= 3 && t_ps <= 5) 
              %           | (t_ps >= 10)
        E_val = 12;
end
E_field(i) = E_val * 1e5; %kV/cm -> V/m
end
figure; %for control only
plot(t*1e12, E_field/1e5, 'LineWidth', 2);
xlabel('Time (ps)');
ylabel('E field (kV/cm)');
title('Applied E-field');
grid on;

for i = 1 : N-1
    %current values
    w_now_j = w(i);
    w_now_eV = w_now_j * 1/q;
    v_now = v(i);
    E_now = E_field(i); %v/m
    
    %bdry check
    if w_now_eV < 0.038
        w_now_eV = 0.038; 
        %upper bdry is discussed in part III.E
    end
    
    tp = calc_tau_p(w_now_eV); %s
    tw = calc_tau_w(w_now_eV); %s
    m_eff = calc_effective_mass(w_now_eV, m0); %kg
    
    %cons. of momentum dv/dt = (qE / m*) - (v / tp)
    dv_dt = (q * E_now / m_eff) - (v_now / tp);

    %cons. of energy dw/dt = qEv - (w - w0)/tw
    dw_dt = (q*E_now * v_now) - ((w_now_j - w(1)) /tw);
    
    %euler
    v(i+1) = v_now + dv_dt * dt;
    w(i+1) = w_now_j + dw_dt * dt;
    
    if w(i+1) < w(1) %checking again, < not possible
       w(i+1) = w(1);
    end
end

%visualization
t_plot = t * 1e12; 
v_plot = v * 100; % m -> cm
w_plot = w / q; %j -> eV

figure;
subplot(2,1,1);
plot(t_plot, v_plot, 'b', 'LineWidth', 2);
xlabel('Time (ps)');
ylabel('Velocity (cm/s)');
title('Electron Velocity');
grid on;
xline([2 3 5 10], 'r'); %e field change pts

subplot(2,1,2);
plot(t_plot, w_plot, 'r', 'LineWidth', 2);
xlabel('Time (ps)');
ylabel('Energy (eV)');
title('Electron Energy');
grid on;
xline([2 3 5 10], 'b');

    

%I wanted to make use of func.s for the calc. of taus
%(valid for 38 meV < w < 350 meV)
function tp = calc_tau_p(w) % w in eV 
    tp_ps = -29058.1 * w^6 + 37288.3 * w^5 ... 
            -18885.1 * w^4 + 4787.52 * w^3 ...
            -635.55  * w^2 + 41.0066 * w - 0.74317;
            
    tp = tp_ps * 1e-12; % ps -> sec
end

function tw = calc_tau_w(w) 
    tw_ps = 5595.75 * w^6 - 9458.1 * w^5 ...
            + 7175.02 * w^4 - 3066.82 * w^3 ...
            + 703.821 * w^2 - 71.1564 * w + 3.80203;
            
    tw = tw_ps * 1e-12; % ps -> sec
end

function m_eff = calc_effective_mass(w, m0) %eV,kg    
    m_ratio = 6.03057 * w^5 - 13.478 * w^4 ...
              + 9.39288 * w^3 - 1.907 * w^2 ...
              + 0.389707 * w + 0.0443495;
              
    m_eff = m_ratio * m0; %in kg
end