function [L,eh] = latencyplot(DATA,groups,varargin)

% function [L, H] = latencyplot(DATA,GROUPS,varargin)
%
% Plots the escape latencies from the given data using the given indices of different groups.
%
% MANDATORY INPUTS:
% -------------------------------------------------------------------------------------------------
%
%   DATA   - DATA is a cell array that contains the raw water-maze data for a given collection of 
%            project files in the study. See 'help readwmdf' for details.
%
%   GROUPS - A cell array of vectors of animal indices. See 'help getgroups' to determine how to
%            build these vectors.
%
% OUTPUT:
% -------------------------------------------------------------------------------------------------
%
%   L - A T x A matrix of latencies, where T is the maximum number of trials and A is the number of
%       animals.
%
%   H - A cell array of handles for the plot.
%
% OPTIONAL INPUTS:
% -------------------------------------------------------------------------------------------------
% Optional inputs to control plotting behaviour can be provided in parameter/value format. The optional inputs are:
%
%   'days' - A cell array of vectors of trial indices to indicate which trials were on different
%            days. Trials on different days are not joined by a line (i.e. a gap appears in the
%            plot). Default = {} (all trials on the same day).
%
%   'line_width' - A scalar value determing the line widths of the plot. Default = 1.
%
%   'marker_size' - A scalar value determing the size of the symbols on the plot. Default = 10.
%
%   'font_size'  - A scalar value determing the size of the font used in the axis. Default = 12.
%
%   'font_name'  - A string defining which font is used in the axis. Default = lucidasans.
%
%   'grp_edgec'  - A cell array determining the colours of the lines for each group. Default =
%                  {'b','r','g','k','m','c','y'}.
%
%   'grp_facec'  - A cell array determining the colours of the symbols for each group. Default =
%                  {'b','r','g','k','m','c','y'}.
%
%   'grp_symbol' - A cell array determining the symbols for each group. Default = {'o','s','^',
%                  '+','x','d','p'}.
%
%   'grp_style'  - A cell array determining the line styles for each group. Default = {'-',':','-.','--'}.
%
%   'tick_dir'   - Whether ticks should be 'in' or 'out. Default = 'out'.
% 
%--------------------------------------------------------------------------------
%
% 02/2013, Frankland Lab (www.franklandlab.com)
%
% Author: Blake Richards
% Contact: blake.richards@utoronto.ca
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARSE THE INPUT

% define the default optional arguments
optargs = struct('days',[],...
                 'line_width',1,...
                 'marker_size',1,...
                 'font_size',12,...
                 'font_name','lucidasans',...
                 'grp_edgec',{{'b','r','g','k','m','c','y'}},...
                 'grp_facec',{{'b','r','g','k','m','c','y'}},...
                 'grp_symbol',{{'o','s','^','+','x','d','p'}},...
                 'grp_style',{{'-',':','-.','--'}},...
                 'tick_dir','out'); 

% get the optional argument names
optnames = fieldnames(optargs);

% get the number of optional arguments
nargs = length(varargin)/2;

% make sure the property/value pairs are input in pairs
if round(length(varargin)/2) ~= nargs
   error('Expecting propertyName/propertyValue pairs after DATA and GROUPS');
end

% step through the optional arguments, check them, and store them
for pair = reshape(varargin,2,[])

	% make it case insensitive
	inpname = lower(pair{1});

	% check whether the name matches a known option
	if any(strmatch(inpname,optnames))
		switch inpname
			case 'days'
				if isa(pair{2},'cell')
					optargs.(inpname) = pair{2};
				else
					error('days must be a cell array');
				end
			case 'line_width'
				if isa(pair{2},'numeric')
					optargs.(inpname) = pair{2};
				else
					error('line_width must be a scalar value');
				end
			case 'marker_size'
				if isa(pair{2},'numeric')
					optargs.(inpname) = pair{2};
				else
					error('marker_size must be a scalar value');
				end
			case 'font_size'
				if isa(pair{2},'numeric')
					optargs.(inpname) = pair{2};
				else
					error('font_size must be a scalar value');
				end
			case 'font_name'
				if isa(pair{2},'char') && ismember(pair{2},listfonts())
					optargs.(inpname) = pair{2};
				else
					error('font_name must be a valid font name (see help listfonts)');
				end
			case 'grp_edgec'
				if isa(pair{2},'cell')
					optargs.(inpname) = pair{2};
				else
					error('grp_edgec must be a cell array');
				end
			case 'grp_facec'
				if isa(pair{2},'cell')
					optargs.(inpname) = pair{2};
				else
					error('grp_facec must be a cell array');
				end
			case 'grp_symbol'
				if isa(pair{2},'cell')
					optargs.(inpname) = pair{2};
				else
					error('grp_symbol must be a cell array');
				end
			case 'grp_style'
				if isa(pair{2},'cell')
					optargs.(inpname) = pair{2};
				else
					error('grp_style must be a cell array');
				end
			case 'tick_dir'
				if isa(pair{2},'char') && ismember(pair{2},{'in','out'})
					optargs.(inpname) = pair{2};
				else
					error('tick_dir must be ''in'' or ''out''');
				end
		end	
	else
		error('%s is not a recognized parameter name',inpname);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE DATA

% get the latencies
L = latency(DATA);

% make the figure
figure();
set(gcf,'Color',[1 1 1],'Position',[300 200 800 400],'PaperPosition',[0 0 8 4]);
axes('Position',[0.1 0.2 0.8 0.7]);
hold on;

% for each day...
if isempty(optargs.days), days = {[1:size(L,1)]}; else, days = optargs.days; end;
for dd = 1:length(days)

  % get the xdata for these trials
	xdata = dd + linspace(0.2,0.8,length(days{dd}));

	% plot the data for each group on this day
	for gg = 1:length(groups)
	
		% calculate the y-values
		ymean = mean(L(days{dd},groups{gg}),2);
		ysem  = std(L(days{dd},groups{gg}),[],2)./sqrt(size(L,2));
		
		% plot the errorbars and format
		eh{gg} = errorbar(xdata,ymean,ysem);
		set(eh{gg},'Color',optargs.grp_facec{mod(gg-1,length(optargs.grp_facec))+1},...
	             'MarkerEdgeColor',optargs.grp_edgec{mod(gg-1,length(optargs.grp_edgec))+1},...
	             'Marker',optargs.grp_symbol{mod(gg-1,length(optargs.grp_symbol))+1},...
	             'LineStyle',optargs.grp_style{mod(gg-1,length(optargs.grp_style))+1},...
	             'MarkerSize',optargs.marker_size,...
	             'LineWidth',optargs.marker_size);

	end
end

% Add labels and tick marks
xlabel('Day','FontName',optargs.font_name,'FontSize',optargs.font_size);
ylabel('Escape latency (s)','FontName',optargs.font_name,'FontSize',optargs.font_size);
set(gca,'XTick',[1:length(days)]+0.5,'XTickLabel',[1:length(days)]);
set(gca,'YTick',[0:10:60],'YTickLabel',[0:10:60]);

% beautify the figure
axis([1 length(days)+1 0 65]);
set(gca,'TickDir',optargs.tick_dir,'LineWidth',optargs.line_width);
set(gca,'FontName',optargs.font_name,'FontSize',optargs.font_size);

