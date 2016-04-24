function plot_ir_v(Median, CIH2, CIL2, shocks, h, T, vert)
periods = 0:h;
%Vec = reshape(reshape((1:T(1,2)^2), T(1,2), T(1,2))', T(1,2)^2, 1)';

Vec = reshape((1:T(1,2)^2), T(1,2), T(1,2));
vec = Vec(:, shocks(1):shocks(end));
if vert == 1
    vec = vec';
    p_fig(1) = T(1,2);
    p_fig(2) = size(shocks, 2);
else
    p_fig(1) = size(shocks, 2);
    p_fig(2) = T(1,2);
end
positions = reshape(vec, T(1,2) * size(shocks, 2),1)';

% normalize the effects:
stand1 = 0.25/Median(14, 1);
stand2 = 0.25/Median(18, 1);
%
% Median(11:15, :) = stand1*Median(11:15, :) ;
% Median(16:20, :) = stand2*Median(16:20, :);
% 
% CIH2(11:15, :) = stand1*CIH2(11:15, :) ;
% CIH2(16:20, :) = stand2*CIH2(16:20, :);
% 
% CIL2(11:15, :) = stand1*CIL2(11:15, :) ;
% CIL2(16:20, :) = stand2*CIL2(16:20, :);

figure('Name','Impulse Responses','NumberTitle','off')
for p = 1 : T(1,2)*size(shocks, 2) %5:T(1,2)^2
                subplot(p_fig(1), p_fig(2), p);
                plot(periods, Median(positions(p),:),'k',  'LineWidth', 1.5);
                hold on;
                plot(periods, CIH2(positions(p),:),'k--',  'LineWidth', 1.5);
                plot(periods, CIL2(positions(p),:),'k--' ,  'LineWidth', 1.5);
                plot(periods, zeros(1, h+1), ':' );
                hold off;
                xlim(gca, [0 h])
end