function PlotImpulseResponses(irf,var_in,info,fig_id)

% PlotImpulseResponses
%
% Usage:
%   PlotImpulseResponses(irf,var_in,info,fig_id);
%   PlotImpulseResponses(irf,var_in,info);
%
% Purpose:
%   Plots impulse responses for a vector autoregression.
%
% Input:
%   irf     n x h+1 x n array, impulse responses
%   info    structure including
%   .conf       1 x c vector, significance levels for confidence bands
%   .area_color 1 x p vector, color of area for last significance level
%   .area_step  double, reduction of color until first significance level
%   .widths     1 x 2 vector, line widths to be used
%   .fsizes     1 x 2 vector, font sizes to be used
%   .names      1 x n vector of strings, variable names
%   var_in  variable inputs;
%    bands      n x h+1 x 2 x c vector, lower and upper asymm. bands
%    sig        n x h+1 vector, IRF standard deviations
%    []         empty matrix
%   fig_id  integer, figure identifier (optional)
%
% Output:
%   none
%
% Author:
%   Markus Kirchner, June 2012

% Get inputs
nvars          = size(irf,1);
horizon        = size(irf,2)-1;
area_step      = info.area_step;
line_width     = info.widths(1);
line_width_alt = info.widths(2);
fsize          = info.fsizes(1);
fsize_alt      = info.fsizes(2);
names          = info.names;

% Confidence or error levels
if isfield(info,'conf')
    levels = info.conf;
elseif isfield(info,'prob')
    levels = info.prob;

end

% Number of error bands
nbands = numel(levels);    
    
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

% Use asymmetric or symmetric bands
if size(var_in,3)>1  
    low_bands = squeeze(var_in(:,:,1,:));
    up_bands  = squeeze(var_in(:,:,2,:));
elseif ~isempty(var_in)
    low_bands = zeros(nvars,horizon+1,nbands);
    up_bands  = zeros(nvars,horizon+1,nbands);
    for i=1:nbands
        p_level          = (1+levels(i))/2;
        band_width       = norminv(p_level,0,1);
        low_bands(:,:,i) = irf - band_width*var_in;
        up_bands(:,:,i)  = irf + band_width*var_in;
    end
    
    low_68 = low_bands(:,:,1);
    up_68  = up_bands(:,:,1);
    
    save('low_band_68','low_68');
    save('up_band_68','up_68');
    
    low_90 = low_bands(:,:,2);
    up_90  = up_bands(:,:,2);
    
    save('low_band_90','low_90');
    save('up_band_90','up_90');
    
    
    
end


        
% Plot impulse responses
if nargin==4
    figure(fig_id);
else
    figure;
end
grey = [0.3,0.3,0.3];
k = 1;
for j=1:nvars
    subplot(ceil(nvars/plot_fac),plot_fac,k)
    if ~isempty(var_in)
        area_color = info.area_color;
        for i=nbands:-1:1
            shadedplot(0:horizon,...
                low_bands(j,:,i),up_bands(j,:,i),...
                area_color);
            area_color = area_color-area_step;
            hold on
        end
    end
    h(1) = plot(0:horizon,zeros(1,horizon+1),'-','Color',grey,...
        'LineWidth',line_width_alt);
    h(2) = plot(0:horizon,irf(j,:),'-b',...
        'LineWidth',line_width);
    title(names(j),'FontSize',fsize)
    if plot_fac==2
        if j==nvars || j==nvars-1
            xlabel('Periods','FontSize',fsize_alt)
        end
    elseif plot_fac==3
        if j==nvars || j==nvars-1 || j==nvars-2
            xlabel('Periods','FontSize',fsize_alt)
        end
    else
        if j==nvars || j==nvars-1 || j==nvars-2 || j==nvars-3
            xlabel('Periods','FontSize',fsize_alt)
        end
    end
    set(gca,'FontSize',fsize_alt)
    xlim([0 horizon])
    box off
    axis tight
    k = k+1;
end
if ~isempty(var_in)
    if length(levels)==1
        legend1 = legend(h(2),['Impulse responses with ',...
            num2str(conf*100) ' percent error bands']);
    elseif length(levels)==2
        legend1 = legend(h(2),['Impulse responses with ',...
            num2str(conf(1)*100) ' and ',...
            num2str(conf(2)*100) ' percent error bands']);
    elseif length(levels)==3
        legend1 = legend(h(2),['Impulse responses with ',...
            num2str(levels(1)*100) ', ',...
            num2str(levels(2)*100) ' and ',...
            num2str(levels(3)*100) ' percent error bands']);
    end
else
    legend1 = legend(h(2),'Point estimates of impulse responses');
end

set(legend1,'FontSize',fsize_alt,...
    'Orientation','horizontal',...
    'Position',[0.327835 0.00079499 0.3672619 0.05176190476],...
    'Box', 'off');