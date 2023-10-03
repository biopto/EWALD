function [sim_time] = show_time(time)

    if time<60
        sim_time = sprintf('%.3f seconds', time);
    elseif time>=60 && time<3600
        sim_time = sprintf('%.0f minutes and %.0f seconds', floor(time/60), round(rem(time, 60)));
    else
        sim_time = sprintf('%.0f hours and %.0f minutes', floor(time/3600), round(rem(time, 3600)/60));
    end
    fprintf('Computation time: %s.\n', sim_time);

end