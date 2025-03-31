function [HO,FS] = process_ISMIP_HOM_ensemble_experiment_D( foldername)
% Read and process the Pattyn et al. (2008) ensemble of ISMIP-HOM
% experiments

ensemble = read_ISMIP_HOM_ensemble_experiment_D( foldername);
model_types = ISMIP_HOM_model_types;

Ls = [5,10,20,40,80,160];

for Li = 1:6

  L = Ls( Li);

  if L<10
    ex = ['L00' num2str(L)];
  elseif L<100
    ex = ['L0'  num2str(L)];
  else
    ex = ['L'   num2str(L)];
  end

  x_FS = linspace(0,1,200)';
  u_FS = [];
  u_HO = [];

  flds = fields(ensemble.(ex));
  for mi = 1:length(flds)
    m = flds{mi};

    x  = ensemble.(ex).(m).x(:,1);
    yi = 1;
    u  = ensemble.(ex).(m).u(:,yi);
 
    % Determine if this model is FS or HO
    is_FS = false;
    is_HO = false;
    for mii = 1:size(model_types,1)
      if strcmpi(model_types{mii,1},m)
        if strcmpi(model_types{mii,2},'FS')
          is_FS = true;
        else
          is_HO = true;
        end
      end
    end
    if ~(is_FS || is_HO)
      for mii = 1:size(model_types,1)
        if strcmpi(model_types{mii,1}(1:3),m(1:3))
          if strcmpi(model_types{mii,2},'FS')
            is_FS = true;
          else
            is_HO = true;
          end
        end
      end
    end
    if ~(is_FS || is_HO)
      % Unknown model?
      continue
    end

    % Add to data ranges for HO/FS models
    up = interp1(x,u,x_FS);
    if is_FS
      u_FS(:,end+1) = up;
    elseif is_HO
      u_HO(:,end+1) = up;
    else
      error('whaa!')
    end

  end

  % Missing data points
  m = true(size(x_FS));
  for i = 1:length(x_FS)
    if sum(isnan(u_FS(i,:)))+sum(isnan(u_HO(i,:)))>0
      m(i) = false;
    end
  end
  u_FS = u_FS(m,:);
  u_HO = u_HO(m,:);
  x_FS = x_FS(m);

  % ISMIP-HOM ensemble data
  FS.(ex).x     = x_FS;
  HO.(ex).x     = x_FS;
  FS.(ex).u_av  = mean(u_FS,2);
  HO.(ex).u_av  = mean(u_HO,2);
  FS.(ex).u_min = u_FS(:,1) + Inf;
  HO.(ex).u_min = u_HO(:,1) + Inf;
  FS.(ex).u_max = u_FS(:,1) - Inf;
  HO.(ex).u_max = u_HO(:,1) - Inf;

  for mi = 1: size( u_FS,2)
    FS.(ex).u_min = min( FS.(ex).u_min, u_FS( :,mi));
    HO.(ex).u_min = min( HO.(ex).u_min, u_HO( :,mi));
    FS.(ex).u_max = max( FS.(ex).u_max, u_FS( :,mi));
    HO.(ex).u_max = max( HO.(ex).u_max, u_HO( :,mi));
  end

end

end