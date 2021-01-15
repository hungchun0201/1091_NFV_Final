format long;
arrival_rate = 4;
service_rate = 5;
capacity = 100;
loss_constraint = 1e-50;
C = 1;

rho = arrival_rate / service_rate;


m_list = 1:capacity;
PK_list = zeros(1, capacity);
cost_list = zeros(1, capacity);
for m = 1:capacity
    sum_from_P0_to_Pm = (1 - rho^m) / (1 - rho);
    sum_from_Pm_to_PK = rho^(m - 1) * rho / 2 * ((1 - (rho / 2)^(capacity - m + 1)) / (1 - rho / 2));
    %fprintf("sum_from_P0_to_Pm: %f\n",sum_from_P0_to_Pm);
    P0 = (sum_from_P0_to_Pm + sum_from_Pm_to_PK)^ - 1;
    %fprintf("P0: %f\n",P0);

    PK = rho^(m - 1) * (rho / 2)^(capacity - m + 1) * P0;
    PK_list(m) = PK;

    cost = (C * sum_from_P0_to_Pm + 2 * C * sum_from_Pm_to_PK) * P0;
    cost_list(m) = cost;
    % fprintf("round: %d, loss rate: %.50f\n",m,PK);
    % fprintf("round: %d, cost: %.10f\n",m,cost);
end

plot_rho_fixed(rho,capacity, m_list, PK_list, cost_list);

function plot_rho_fixed(rho,capacity, m_list, PK_list, cost_list)
    clf

    % ---------------------------------------------------------------------------- %
    %                              Determine Optimal m                             %
    % ---------------------------------------------------------------------------- %

    derivate_cost = diff(cost_list) ./ diff(m_list);
    x_derivative = m_list(2:end);
    optimal_m = 2;
    for i = 2:m_list(end)
        if (derivate_cost(i) >- 0.0005)
            optimal_m = i;
            break
        end
    end

    % ---------------------------------------------------------------------------- %
    %                                  Plot Figure                                 %
    % ---------------------------------------------------------------------------- %

    subplot(1, 2, 1);
    semilogy(m_list, PK_list);

    hold on
    % xline(optimal_m);
    xlabel("m");
    ylabel("Packet Loss Probability")
    title(['Packet Loss Probability v.s. m, ρ(λ/μ) = ', num2str(rho),', capacity = ',num2str(capacity)]);
    hold on

    % dim = [.2 .5 .3 .3];
    % str = ['The packet loss prob is ', num2str(PK_list(optimal_m)), ', if m = ', num2str(optimal_m)];
    % annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');

    subplot(1, 2, 2);

    plot(m_list, cost_list);
    hold on
    plot(x_derivative, derivate_cost);
    hold on
    % xline(optimal_m);

    % dim = [.6 .5 .3 .3];
    % str = ['The cost is ', num2str(cost_list(optimal_m)), ', if m = ', num2str(optimal_m)];
    % annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');

    xlabel("m");
    ylabel("Cost")
    legend("C(m)", "differential of C(m)", ['Optimal m = ', num2str(optimal_m)])
    title(['Cost function v.s. m, ρ(λ/μ) = ', num2str(rho),', capacity = ',num2str(capacity)]);
end
