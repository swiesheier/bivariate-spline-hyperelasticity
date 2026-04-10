function data = get_data(type, model, modes, varargin)
% GET_DATA Load data for a given type, model, and modes.

switch lower(type)

    case 'treloar'
        data = struct();
        data.type  = "treloar";
        data.modes = upper(string(modes));
        for i = 1:numel(modes)
            mode = upper(modes{i});
            data.(mode) = load_treloar(mode);
        end

    otherwise
        error('Unknown data type');
end
end
