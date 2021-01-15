format long;
arrival_rate = 4;
service_rate = 5;
capacity = 2000;
loss_constraint = 1e-50;
C = 1;

rho = arrival_rate/service_rate;

for capacity = [100,500,1000,5000]
    rho_list = 0.101:0.1:1.901;
    optimal_m_list = zeros(1,length(rho_list));
    loss_prob_list = zeros(1,length(rho_list));
    a = 0;
    for rho = rho_list
        a = a+1;
        m_list = 1:capacity;
        PK_list = zeros(1,capacity);
        cost_list = zeros(1,capacity);
        for m=1:capacity
            sum_from_P0_to_Pm = (1-rho^m)/(1-rho);
            sum_from_Pm_to_PK = rho^(m-1)*rho/2*((1-(rho/2)^(capacity-m+1))/(1-rho/2));
            %fprintf("sum_from_P0_to_Pm: %f\n",sum_from_P0_to_Pm);
            P0 = (sum_from_P0_to_Pm + sum_from_Pm_to_PK)^-1;
            %fprintf("P0: %f\n",P0);
                
            PK = rho^(m-1)*(rho/2)^(capacity-m+1)*P0;
            PK_list(m) = PK;
            
            cost = (C*sum_from_P0_to_Pm+2*C*sum_from_Pm_to_PK)*P0;
            cost_list(m) = cost;
            % fprintf("round: %d, loss rate: %.50f\n",m,PK);
            % fprintf("round: %d, cost: %.10f\n",m,cost);
        end

        derivate_cost = diff(cost_list)./diff(m_list);
        x_derivative = m_list(2:end);
        optimal_m = 2;
        for i=2:m_list(end)
            if(derivate_cost(i-1)>-0.0005)
                optimal_m = i;
                break
            end
        end
        
        
        optimal_m_list(a) = optimal_m;
        loss_prob_list(a) = PK_list(optimal_m);

        %plot_rho_fixed(rho,m_list,PK_list,cost_list);
        
    end
    subplot(1,2,1);
    plot(rho_list,optimal_m_list);
    hold on
    xlabel("ρ");
    ylabel("Optimal m");
    title("The determined optimal m for different ρ and capacity");


    subplot(1,2,2);
    semilogy(rho_list,loss_prob_list);
    xlabel("ρ");
    ylabel("Packet Loss Probability under the optimal m");
    title("The packet loss probability under optimal m for different ρ and capacity");
    hold on
end
subplot(1,2,1);
legend("Capacity = 100","Capacity = 500","Capacity = 1000","Capacity = 5000");
subplot(1,2,2);
legend("Capacity = 100","Capacity = 500","Capacity = 1000","Capacity = 5000");
function plot_rho_fixed(rho,m_list,PK_list,cost_list)
    clf

% ---------------------------------------------------------------------------- %
%                              Determine Optimal m                             %
% ---------------------------------------------------------------------------- %

    derivate_cost = diff(cost_list)./diff(m_list);
    x_derivative = m_list(2:end);
    optimal_m = 2;
    for i=2:m_list(end)
        if(derivate_cost(i)>-0.0005)
            optimal_m = i;
            break
        end
    end

% ---------------------------------------------------------------------------- %
%                                  Plot Figure                                 %
% ---------------------------------------------------------------------------- %

    subplot(1,2,1);
    semilogy(m_list,PK_list);
    
    hold on
    xline(optimal_m);
    xlabel("m");
    ylabel("Packet Loss Probability")
    title(['Packet Loss Probability v.s. m, ρ(λ/μ) = ',num2str(rho) ]);
    hold on
    %dy=diff(PK_list)./diff(m_list)
    %plot(m_list(2:end),dy)
    dim = [.2 .5 .3 .3];
    str = ['The packet loss prob is ',num2str(PK_list(optimal_m)),', if m = ',num2str(optimal_m)];
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
    subplot(1,2,2);

    plot(m_list,cost_list);
    hold on
    plot(x_derivative,derivate_cost);
    hold on
    xline(optimal_m);
    
    dim = [.6 .5 .3 .3];
    str = ['The cost is ',num2str(cost_list(optimal_m)),', if m = ',num2str(optimal_m)];
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
    xlabel("m");
    ylabel("Cost")
    legend("C(m)","derivate of C(m)",['Optimal m = ',num2str(optimal_m)])
    title(['Cost function v.s. m, ρ(λ/μ) = ', num2str(rho)]);
end