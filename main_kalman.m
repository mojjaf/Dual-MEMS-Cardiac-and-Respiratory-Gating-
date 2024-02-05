function [] = main_kalman(amp,tmp,fs)
    % time step
    
    samples=1:numel(amp);
    t=(samples/fs)';
    dt = 1/fs;

%     t=(0:dt:2)';
    n = numel(t);

    %ground truth
    signal = amp; %%% RPM Signal

    % state matrix
    X = zeros(2,1);

    % covariance matrix
    P = zeros(2,2);

    % kalman filter output through the whole time
    X_arr = zeros(n, 2);

    % system noise
    Q = [0.04 0;
         0 1];

    % transition matrix
    F = [1 dt;
         0 1]; 

    % observation matrix 
    H = [1 0];
    
    s1 = tmp-mean(tmp);
%     s1 = zscore(tmp);

    if length(signal)>length(s1)
        padsize=length(signal)-length(s1);
        s1=padarray(s1',[0,padsize],0,'post');
    else
        padsize=length(s1)-length(signal);
        signal=padarray(signal,[0,padsize],0,'post');
    end

    % variance of signal 1 
%     s1_var = 0.08*ones(size(t)); 
%     s1 = generate_signal(signal, s1_var);
    s1_var = var(s1).*ones(length(s1),1); 
    s1=[s1', s1_var];
    s1(:,1)=s1(:,1)+signal;

    % variance of signal 2 
%     s2_var = 0.01*(cos(8*t)+10*t);
%     s2 = generate_signal(signal, s2_var);
% 
%     % variance of signal 3 
%     s3_var = 0.02*(sin(2*t)+2);
%     s3 = generate_signal(signal, s3_var);
% 
%     % variance of signal 4 
%     s4_var = 0.06*ones(size(t)); 
%     s4 = generate_signal(signal, s4_var);

    % fusion
    for i = 1:n
        if (i == 1)
            [X, P] = init_kalman(X, s1(i, 1)); % initialize the state using the 1st sensor
        else
            [X, P] = prediction(X, P, Q, F);

            [X, P] = update(X, P, s1(i, 1), s1(i, 2), H);
%             [X, P] = update(X, P, s2(i, 1), s2(i, 2), H);
%             [X, P] = update(X, P, s3(i, 1), s3(i, 2), H);
%             [X, P] = update(X, P, s4(i, 1), s4(i, 2), H);
        end

        X_arr(i, :) = X;
    end
    figure
    plot(t, signal, 'LineWidth', 4);
    hold on;
    plot(t, s1(:, 1), '--', 'LineWidth', 1);
%     plot(t, s2(:, 1), '--', 'LineWidth', 1);
%     plot(t, s3(:, 1), '--', 'LineWidth', 1);
%     plot(t, s4(:, 1), '--', 'LineWidth', 1);
    plot(t, X_arr(:, 1), 'LineWidth', 4);
    hold off;
    grid on;
    legend('Ground Truth', 'Sensor Input 1', 'Fused Output');
end

function [s] = generate_signal(signal, var)
    noise = randn(size(signal)).*sqrt(var);

    s(:, 1) = signal + noise;
    s(:, 2) = var; 
end

function [X, P] = init_kalman(X, y)
    X(1,1) = y;
    X(2,1) = 0;

    P = [100 0;
         0   300];
end

function [X, P] = prediction(X, P, Q, F)
    X = F*X;
    P = F*P*F' + Q;
end

function [X, P] = update(X, P, y, R, H)
    Inn = y - H*X;
    S = H*P*H' + R;
    K = P*H'/S;

    X = X + K*Inn;
    P = P - K*H*P;
end