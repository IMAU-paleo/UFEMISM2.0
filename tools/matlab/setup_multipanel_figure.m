function H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver)
% Set up a single figure with mulitple axeses

%% Safety

% margins_hor and margins_ver need at least two values each
if length( margins_hor) < 2
  error('need two values (left and right) for margins_hor!')
end
if length( margins_ver) < 2
  error('need two values (top and bottom) for margins_ver!')
end

% If a single value for wa or ha is provided, make all axeses the same size
if length( wa) == 1
  wa = zeros( 1, length( margins_hor) - 1) + wa;
else
  if length( wa) ~= length( margins_hor) - 1
    error('length( wa) should be length( margins_hor) - 1!')
  end
end
if length( ha) == 1
  ha = zeros( 1, length( margins_ver) - 1) + ha;
else
  if length( ha) ~= length( margins_ver) - 1
    error('length( ha) should be length( margins_ver) - 1!')
  end
end

%% Set up GUI

% Number of axeses
nax = length( margins_hor) - 1;
nay = length( margins_ver) - 1;

% Figure size
wf = sum( margins_hor) + sum( wa);
hf = sum( margins_ver) + sum( ha);

% Create figure
H.Fig = figure('position',[20,20,wf,hf],'color','w');

% Create axeses
for i = 1: nay
  for j = 1: nax
    x = sum( margins_hor( 1:j)) + sum( wa( 1:j-1));
    y = sum( margins_ver( 1:i)) + sum( ha( 1:i-1));
    yb = hf - y - ha( i);
    H.Ax{ i,j} = axes('parent',H.Fig,'units','pixels','position',[x,yb,wa( j),ha( i)],'fontsize',20,...
      'xgrid','on','ygrid','on');
  end
end

end