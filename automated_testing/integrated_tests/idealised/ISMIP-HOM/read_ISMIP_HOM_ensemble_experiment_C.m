function ensemble = read_ISMIP_HOM_ensemble_experiment_C( foldername)
% Read Pattyn et al. (2008) ensemble results

models = dir(foldername);
models = models(3:end);

ensemble.L160 = [];
ensemble.L080 = [];
ensemble.L040 = [];
ensemble.L020 = [];
ensemble.L010 = [];
ensemble.L005 = [];

for mi = 1: length(models)

  modeldata = dir([foldername '/' models(mi).name]);
  modeldata = modeldata(3:end);

  % Go over all experiments, check if this model has them.
  flds = fields(ensemble);
  for xi = 1:length(flds)
    ex = flds{xi};

    for di = 1:length(modeldata)
      mdname = modeldata(di).name;
      str = [strrep(ex,'L','c') '.txt'];
      if length(mdname) >= length(str)
        if strcmpi(mdname(end-length(str)+1:end),str)
          % This is the experiment from this model

          disp(['Reading data from model ' models(mi).name ', experiment ' strrep(ex,'L','c')])

          fid = fopen([foldername '/' models(mi).name '/' mdname]);
          temp = textscan(fid,'%s','delimiter','\n','MultipleDelimsAsOne',1); temp = temp{1};
          fclose(fid);

          n = length(temp);
          nx = sqrt(n);
          if nx-floor(nx)>0
            error('whaa!')
          end
          x_vec = zeros(n,1);
          y_vec = zeros(n,1);
          u_vec = zeros(n,1);
          v_vec = zeros(n,1);
          w_vec = zeros(n,1);
          for i = 1:n
            temp2 = textscan(temp{i},'%f %f %f %f %f %f %f %f %f %f');
            x_vec(i) = temp2{1};
            y_vec(i) = temp2{2};
            u_vec(i) = temp2{3};
            v_vec(i) = temp2{4};
            w_vec(i) = temp2{5};
          end

          ensemble.(ex).(models(mi).name).x = reshape(x_vec,[nx,nx]);
          ensemble.(ex).(models(mi).name).y = reshape(y_vec,[nx,nx]);
          ensemble.(ex).(models(mi).name).u = reshape(u_vec,[nx,nx]);
          ensemble.(ex).(models(mi).name).v = reshape(v_vec,[nx,nx]);
          ensemble.(ex).(models(mi).name).w = reshape(w_vec,[nx,nx]);

          if (ensemble.(ex).(models(mi).name).x(1,1) == ensemble.(ex).(models(mi).name).x(end,1))
            ensemble.(ex).(models(mi).name).x = ensemble.(ex).(models(mi).name).x';
            ensemble.(ex).(models(mi).name).y = ensemble.(ex).(models(mi).name).y';
            ensemble.(ex).(models(mi).name).u = ensemble.(ex).(models(mi).name).u';
            ensemble.(ex).(models(mi).name).v = ensemble.(ex).(models(mi).name).v';
            ensemble.(ex).(models(mi).name).w = ensemble.(ex).(models(mi).name).w';
          end

        end
      end
    end
  end

end

end