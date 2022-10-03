function PlotRegressionResults(data,results,info,fig_id1,fig_id2)

% PlotRegressionResults
%
% Usage:
%   PlotRegressionResults(data,results,info,fig_id1,fig_id2);
%   PlotRegressionResults(data,results,info,fig_id1);
%   PlotRegressionResults(data,results,info);
%
% Purpose:
%   Plots data against fitted values as well as residuals from a vector 
%   autoregression, including one standard error bands for residuals.
%
% Input:
%   data    T+p x n matrix, data with obs./variables in rows/columns
%   results structure including
%   .fitted     T x n matrix, fitted values
%   .resid      T x n matrix, regression residuals
%   .Omega      n x n matrix, residual variance-covariance matrix
%   info    structure including
%   .lag_order  integer, VAR lag order p
%   .dates      T+p x 1 vector, dates corresponding to data set
%   .names      1 x n vector of strings, variable names
%   .widths     1 x 2 vector, line widths to be used
%   .fsizes     1 x 2 vector, font sizes to be used
%   fig_id1 integer, figure identifier data vs. fitted values (optional)
%   fig_id2 integer, figure identifier residuals (optional)
%
% Output:
%   none
%
% Author:
%   Markus Kirchner, May 2012

% Get inputs
[nobs nvars]   = size(data);
line_width     = info.widths(1);
line_width_alt = info.widths(2);
fsize          = info.fsizes(1);
fsize_alt      = info.fsizes(2);
names          = info.names;
lag_order      = info.lag_order;
dates          = info.dates;
fitted         = results.fitted;
resid          = results.resid;
std_error      = sqrt(diag(results.Omega));

% Prepare plots
if nvars<=3
    plot_fac = 1;
elseif nvars>3 && nvars<=8
    plot_fac = 2;
elseif nvars>8 && nvars<13
    plot_fac = 3;
else
    plot_fac = 4;
end
        
% Plot data against fitted values
if nargin>3
    figure(fig_id1);
else
    figure;
end
k = 1;
for j=1:nvars
    subplot(ceil(nvars/plot_fac),plot_fac,k)
    plot(dates,data(:,j),'-b',...
        'LineWidth',line_width);
    hold on
    plot(dates,[NaN(lag_order,1); fitted(:,j)],'--r',...
        'LineWidth',line_width);
    title(names(j),'FontSize',fsize)
    set(gca,'FontSize',fsize_alt)
    xlim([dates(1) dates(end)])
    axis tight
    k = k+1;
end
legend1 = legend('Data','Fitted values');
set(legend1,'FontSize',fsize_alt,...
    'Orientation','horizontal',...
    'Position',[0.327835 0.00079499 0.3672619 0.05176190476],...
    'Box', 'off');

% Plot regression residuals
if nargin>4
    figure(fig_id2);
else
    figure;
end
grey = [0.3,0.3,0.3];
k = 1;
for j=1:nvars
    subplot(ceil(nvars/plot_fac),plot_fac,k)
    h(1) = plot(dates,zeros(1,nobs),'-','Color',grey,...
        'LineWidth',line_width_alt);
    hold on
    h(2) = plot(dates,[NaN(lag_order,1); resid(:,j)],'-b',...
        'LineWidth',line_width);
    h(3) = plot(dates,std_error(j)*ones(1,nobs),'--','Color',grey,...
        'LineWidth',line_width_alt);
    h(4) = plot(dates,-std_error(j)*ones(1,nobs),'--','Color',grey,...
        'LineWidth',line_width_alt);
    title(names(j),'FontSize',fsize)
    set(gca,'FontSize',fsize_alt)
    xlim([dates(1) dates(end)])
    axis tight
    k = k+1;
end
legend1 = legend(h([2 3]),'Residuals','+/- One standard error');
set(legend1,'FontSize',fsize_alt,...
    'Orientation','horizontal',...
    'Position',[0.327835 0.00079499 0.3672619 0.05176190476],...
    'Box', 'off');