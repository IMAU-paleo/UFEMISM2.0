function analyse_integrated_test( varargin)
  if ~isempty(varargin)
    analyse_integrated_test_H_dHdt_flowline(     varargin{1});
    analyse_integrated_test_H_dHdt_local(        varargin{1});
    analyse_integrated_test_H_u_flowline(        varargin{1});
    analyse_integrated_test_dHdt_invfric_invBMB( varargin{1});
  else
    analyse_integrated_test_H_dHdt_flowline;
    analyse_integrated_test_H_dHdt_local;
    analyse_integrated_test_H_u_flowline;
    analyse_integrated_test_dHdt_invfric_invBMB;
  end
end