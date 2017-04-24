function [ExpParam, DGZ] = LoadDGZdata(DGZfile)

%============================= LoadDGZdata.m ==============================
% Read read data for a single experimental block from the input .dgz file.

  
DGZ                     = dg_read(DGZfile);                                 % Read data
par_num                 = size(DGZ.e_pre,1);                            	% Count how many 
ExpParam.DGZfilename    = DGZfile;                                          % Get experiment filename
ExpParam.Experiment     = DGZ.e_pre{1}{2};                                  % Get experiment name
ExpParam.Subject        = DGZ.e_pre{2}{2};                                  % Get subject name
try
    ExpParam.ADCdegHV    = DGZ.e_params{1}{2};                              % Get ADC/deg scale
    for j = 3:2:(floor(par_num/2)*2)                                    	% For each parameter...
        DGZ.e_pre{j}{2}(ismember(DGZ.e_pre{j}{2}, ' ()!?*-:%')) = '_';   	% Remove punctuation from parameter names
        if ~isempty(DGZ.e_pre{j}{2})                                       	% If parameter name cell is not empty...
            try
                eval(sprintf('ExpParam.%s = %d;', DGZ.e_pre{j}{2}, str2double(DGZ.e_pre{j+1}{2})));
            end
        else
            break;                                                              
        end
    end                                       
%     Eye(i).QNXPosX = DGZ.ems{1}{2};                                     	% Horizontal eye signal
%     Eye(i).QNXPosY = DGZ.ems{1}{3};                                     	% Vertical eye signal
end
