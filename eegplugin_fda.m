% eegplugin_fda() - Plugin to perform Functional Data Analysis functions
%                   called by EEGLAB at startup to create a menu.
%
% Usage:
%   >> eegplugin_fda(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Authors: Mohammad Fayaz
%
% This structure is from eegplugin_erpsource.m 
function vers = eegplugin_fda(fig, trystrs, catchstrs)

vers = 'FDA1.0'; 
if nargin < 3
    error('eegplugin_roiconnect requires 3 arguments');
end

toolsmenu = findobj(fig, 'tag', 'tools');

% we create the menu below
cb2 = [ 'try, LASTCOM = pop_fda_2(EEG);' catchstrs.add_to_hist  ]; 
cb3 = [ 'try, LASTCOM = pop_fda_3(EEG);' catchstrs.add_to_hist  ]; 
cb4 = [ 'try, LASTCOM = pop_fda_4(EEG);' catchstrs.add_to_hist  ]; 
%roi_m = uimenu( toolsmenu, 'label', 'FDA', 'CallBack', cb1, 'separator', 'on');
roi_m0 = uimenu( toolsmenu, 'label', 'FDA');
roi_m00 = uimenu( roi_m0, 'label', 'Epoch Anlaysis');
roi_m01 = uimenu( roi_m0, 'label', 'ERP Analysis');
roi_m00_2 = uimenu( roi_m00, 'label', 'FPCA', 'CallBack', cb2, 'separator', 'off');
roi_m01_1 = uimenu( roi_m01, 'label', 'Descriptive', 'CallBack', cb4, 'separator', 'off');
roi_m01_2 = uimenu( roi_m01, 'label', 'FCCA', 'CallBack', cb3, 'separator', 'off');

end