function ensemble = read_ISMIP_HOM_ensemble_experiment_B( foldername)
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
      str = [strrep(ex,'L','b') '.txt'];
      if length(mdname) >= length(str)
        if strcmpi(mdname(end-length(str)+1:end),str)
          % This is the experiment from this model

          disp(['Reading data from model ' models(mi).name ', experiment ' strrep(ex,'L','b')])

          fid = fopen([foldername '/' models(mi).name '/' mdname]);
          temp = textscan(fid,'%s','delimiter','\n','MultipleDelimsAsOne',1); temp = temp{1};
          fclose(fid);

          n = length(temp);
          x_vec = zeros(n,1);
          u_vec = zeros(n,1);
          w_vec = zeros(n,1);

          for i = 1:n
            temp2 = textscan(temp{i},'%f %f %f %f %f');
            x_vec(i) = temp2{1};
            u_vec(i) = temp2{2};
            w_vec(i) = temp2{3};
          end

          ensemble.(ex).(models(mi).name).x = x_vec;
          ensemble.(ex).(models(mi).name).u = u_vec;
          ensemble.(ex).(models(mi).name).w = w_vec;

        end
      end
    end
  end

end

end