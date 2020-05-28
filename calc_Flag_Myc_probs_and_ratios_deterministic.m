function [Flag_probs, Flag_probs_given_nucl, FM_dyn] = calc_Flag_Myc_probs_and_ratios_deterministic(M_without_probs, ts, Flag_assembly_prob_fun, ...
    myc_indeces, flag_indeces, init_FM_state, FM_states)
% solves the ODEs to model the Flag/Myc histone dynamics experiments

    my_ode_options = odeset('NonNegative',1,'MaxStep',0.01);

    Flag_Myc_odefun_temp = @(t,y) Flag_Myc_odefun(t, y, M_without_probs, Flag_assembly_prob_fun, myc_indeces, flag_indeces);

    [ts_temp, FM_dyn] = ode23(Flag_Myc_odefun_temp,ts,init_FM_state,my_ode_options);  % ode45 may be more accurate, but tests showed no difference and ode23 is a bit faster
    if any(FM_dyn<-1e-3)
        error("concerning negative FM_state probs found in ODE solution")
    end
    if any(FM_dyn<0)
        warning("negative FM_state probs found in ODE solution")
        FM_dyn(FM_dyn<0) = 0;
    end
    if ~all(ts_temp == ts')
        error("ode ts do not agree with input ts")
    end
    if any(isnan(FM_dyn))
        error("ode solution contains NaN")
    end

    Flag_probs = [sum(FM_dyn(:,FM_states(:,1) == 2),2) sum(FM_dyn(:,FM_states(:,2) == 2),2) sum(FM_dyn(:,FM_states(:,3) == 2),2)];
    Flag_probs_given_nucl = Flag_probs ./ [sum(FM_dyn(:,FM_states(:,1) > 0),2) sum(FM_dyn(:,FM_states(:,2) > 0),2) sum(FM_dyn(:,FM_states(:,3) > 0),2)];
    
%     T = zeros(8,27);  % to transform Flag/Myc configuration distributions into 'any nucleosome' configuration distributions
%     T(1,1:8) = 1;
%     T(2,9:12) = 1;
%     T(3,13:16) = 1;
%     T(4,17:20) = 1;
%     T(5,21:22) = 1;
%     T(6,23:24) = 1;
%     T(7,25:26) = 1;
%     T(8,27) = 1;
%     
%     figure
%     subplot(5,1,1)
%     plot(ts_temp, Flag_assembly_prob_fun(ts_temp))
% 
%     subplot(5,1,2)
%     plot(ts_temp, FM_dyn)
% 
%     subplot(5,1,3)
%     plot(ts_temp, T * FM_dyn')
% 
%     subplot(5,1,4)
%     plot(ts_temp, Flag_probs)
% 
%     subplot(5,1,5)
%     plot(ts_temp, Flag_probs_given_nucl)
end
