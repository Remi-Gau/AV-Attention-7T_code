function clim = set_clim(condition)
% returns some pre-specified values for clim.max and clim.min for plotting laclim.minar
% profiles and s-parameters insets
%
% clim.min.profile : min for laminar profiles
% clim.max.profile : max for laminar profiles
% clim.min.inset : min for s-parameters insets ; first column for constant and the second for linear
% clim.max.inset : max for s-parameters insets ; first column for constant and the second for linear
%
% condition
% 0 - (default) no pre-specified limit: use the data from the graph to specify limits
% 1 - for activations
% 2 - for deactivations in A1 and PT
% 3 - for deactivations in V1-2-3
% 4 - for cross modal effects for A1-PT
% 5 - for cross modal effects for V1-2-3
% 6 - for attention effects
% 11 - for MVPA

if nargin<1 || isempty(condition)
    condition = 0;
end

switch condition
    case 0
        clim = [];
    case 1
        % % for activations
        clim.min.profile = -.1;
        clim.max.profile = 3;
        clim.min.inset = [-.1 -.02];
        clim.max.inset = [2.7 .78];
    case 2
        % for deactivations in A1 and PT
        clim.min.profile = -.3;
        clim.max.profile = .1;
        clim.min.inset = [-.5 -.15];
        clim.max.inset = [.5 .25];
    case 3
        % for deactivations in V1-2-3
        clim.min.profile = -1.6;
        clim.max.profile = .1;
        clim.min.inset = [-2 -.5];
        clim.max.inset = [.5 .1];
    case 4
        % for cross modal effects for A1-PT
        clim.min.profile = -.22;
        clim.max.profile = .6;
        clim.min.inset = [-1.2 -.3];
        clim.max.inset = [1.2 .3];
    case 5
        % for cross modal effects for V1-2-3
        clim.min.profile = -.31;
        clim.max.profile = .1;
        clim.min.inset = [-1.2 -.3];
        clim.max.inset = [1.2 .3];
    case 6
        % for attention effects
        clim.min.profile = -.2;
        clim.max.profile = .5;
        clim.min.inset = [-2 -.5];
        clim.max.inset = [.5 .1];
        
    case 11
        % for MVPA
        clim.min.profile = .45;
        clim.max.profile = .77;
        clim.min.inset = [-.2 -.08];
        clim.max.inset = [.4 .08];

end

end

