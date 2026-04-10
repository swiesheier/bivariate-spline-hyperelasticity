function D = load_treloar(mode)
% mode = 'UT' | 'PS' | 'ET'
mode = upper(mode);

fname = fullfile('data','experimental','treloar', sprintf('treloar_%s.csv', mode));
T = readtable(fname);

D = struct();
D.lambda = T.lambda;
D.P11    = T.P_MPa;
D.dt     = [];   % no time data for Treloar
D.w      = 1.0;  % optional weight
end
