% pop_fda_3 - perform ERP analysis with functional canonical correlation   
%             analysis (FCCA). It also smooth curves with B-Spline and 
%             Fourier basis functions.The estimation of penalty term is 
%             obtained with generalized cross validavtion (GCV). 
%             It also compute and plot cross-covariance and
%             cross-correlation between two channels. 
%
% Usage:
%  pop_fda_3(EEG); % pop up window asking users to select method
%
% Inputs:
%  EEG - EEGLAB dataset where ICA weights are estimated
%        and they are epoched. 
%
% Output:
%  In version 0.1, this function does not return any output expect a string to 
%  be added to EEGLAB history. 
%
% Author: Mohammad Fayaz. The FDA plug-in GUI codes.
%         In version 0.1, we call functions with permission from fdaM
%                         It is from J.O. Ramsay,Giles Hooker and Spencer Graves, 
%                                    “Functioanl Data Analysis with R and MATLAB”, 
%                                    Springer, 2009 
%                                    Website: : https://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/Matlab/fdaM.zip
%                                    The alternative of this package is in R CRAN,
%                                    called fda (URL:https://cran.r-project.org/web/packages/fda/index.html)
%                                                License:	GPL-2 | GPL-3 [expanded from: GPL (? 2)]                         
%
% Example:
%   % Examples are already useful for users
%   pop_fda_3(EEG); % pop up window asking users to select method
%
% Copyright (C) 2021 Mohammad Fayaz
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.



function com = pop_fda_3(EEG, varargin)
%%%% Tools -> ERP Analysis -> FCCA
fprintf('ERP Analysis \n');
fprintf('Functional Canonical Correlation Analysis (FCCA) \n');


com = '';
if nargin < 1
    help pop_sourcereconstruction;
    return
end

if nargin < 2
    %%%% GUI 
    g = [1 0.5 0.5];
    geometry = { g g g g g g g g g g g g g g g g g g g g g g g g g};
    uilist = { ...
      { 'Style', 'text', 'string', 'Channel Selection 1', 'fontweight', 'bold' } {} {} ...                       % 1
      { 'Style', 'text', 'string', '   Channel number', 'FontAngle', 'italic'  } ...                             % 2
      { 'Style', 'edit', 'string', '1' 'tag' 'SelectedChanIndex_1'} {} ...
      ...
      { 'Style', 'text', 'string', '   Time limits [min max] (msec)', 'FontAngle', 'italic' } ...                % 3
      { 'Style', 'edit', 'string', '-0.2 0.800 ' 'tag' 'TIS_1' } {} ...
      ...
      { 'Style', 'text', 'string', 'Channel Selection 2', 'fontweight', 'bold' } {} {} ...                       % 4
      { 'Style', 'text', 'string', '   Channel number', 'FontAngle', 'italic'  } ...                             % 5
      { 'Style', 'edit', 'string', '5' 'tag' 'SelectedChanIndex_2'} {} ...
      ...
      { 'Style', 'text', 'string', '   Time limits [min max] (msec)', 'FontAngle', 'italic' } ...                % 6         
      { 'Style', 'edit', 'string', '-0.2 0.800 ' 'tag' 'TIS_2' } {} ...
      ...    
      { 'Style', 'text', 'string', 'Preprocessing (Channel 1) ', 'fontweight', 'bold' } {} {} ...                % 7               
      { 'Style', 'text', 'string', '   Choose Basis Function (Default: B-Spline)', 'FontAngle', 'italic' } ...   % 8
      { 'Style', 'popupmenu', 'string', 'B-Spline|Fourier' , 'tag' 'TypeBasisFunc_1' } {} ...
      { 'Style', 'text', 'string', '   Number of Basis', 'FontAngle', 'italic' } ...                             % 9
      { 'Style', 'edit', 'string', '120' 'tag' 'NB_1'  } {} ...
      ...
      { 'Style', 'text', 'string', '   Order of Basis','FontAngle', 'italic' } ...                               % 10
      { 'Style', 'edit', 'string', '6' 'tag' 'NORD_1' } {} ...
      { 'Style', 'text', 'string', '   Parameter Estimation', 'FontAngle', 'italic' } ...                        % 11
      { 'Style', 'checkbox', 'string' 'GCV' 'value' 1 'tag' 'TickGCV_1' } ...
      { 'Style', 'checkbox', 'string' 'Plot GCV' 'value' 1 'tag' 'PlotGCV_1' } ...
      { 'Style', 'text', 'string', 'Preprocessing (Channel 2) ', 'fontweight', 'bold' } {} {} ...                % 12
      { 'Style', 'text', 'string', '   Choose Basis Function (Default: B-Spline)', 'FontAngle', 'italic' } ...   % 13
      { 'Style', 'popupmenu', 'string', 'B-Spline|Fourier' , 'tag' 'TypeBasisFunc_2' } {} ...
      { 'Style', 'text', 'string', '   Number of Basis', 'FontAngle', 'italic' } ...                             % 14
      { 'Style', 'edit', 'string', '120' 'tag' 'NB_2'  } {} ...
      ...
      { 'Style', 'text', 'string', '   Order of Basis','FontAngle', 'italic' } ...                               % 15
      { 'Style', 'edit', 'string', '6' 'tag' 'NORD_2' } {} ...
      { 'Style', 'text', 'string', '   Parameter Estimation', 'FontAngle', 'italic' } ...                        % 16
      { 'Style', 'checkbox', 'string' 'GCV' 'value' 1 'tag' 'TickGCV_2' } ...
      { 'Style', 'checkbox', 'string' 'Plot GCV' 'value' 1 'tag' 'PlotGCV_2' } ...
      { 'Style', 'text', 'string', 'Functional Canonical Correlation Analysis (FCCA)', 'fontweight', 'bold' } {} {} ...      % 17
      { 'Style', 'text', 'string', '   FCCA Status',  'FontAngle', 'italic' } ...                                            % 18
      { 'Style', 'checkbox', 'string' 'Done' 'value' 1 'tag' 'FCCA_STATUS' } {} ...
      { 'Style', 'text', 'string', '   Number of FCCA', 'FontAngle', 'italic' } ...                                          % 19
      { 'Style', 'edit', 'string', '4' 'tag' 'FCCA_NCAN'} {} ...
      { 'Style', 'text', 'string' '   Plot Pairwaise Canonical Weight Function' 'FontAngle', 'italic' } ...               % 20
      { 'Style', 'edit', 'string', '4' 'tag' 'FCCA_NCAN_P'} {} ...
      { 'Style', 'text', 'string', '   Plot Weight Functions', 'FontAngle', 'italic' } ...
      { 'Style', 'checkbox', 'string' 'Diagonal' 'value' 1 'tag' 'FCCA_Diagonal' }  ...
      { 'Style', 'checkbox', 'string' 'All' 'value' 0 'tag' 'FCCA_offDiagonal' }  ...
      { 'Style', 'text', 'string', '   Plot Scores of Each Trials', 'FontAngle', 'italic' } ...
      { 'Style', 'checkbox', 'string' 'Diagonal' 'value' 1 'tag' 'FCCA_Diagonal_Scores' }  ...
      { 'Style', 'checkbox', 'string' 'All' 'value' 0 'tag' 'FCCA_offDiagonal_Scores' }  ...
      { 'Style', 'text', 'string', '   Plot Correlations', 'FontAngle', 'italic' } ...
      { 'Style', 'checkbox', 'string' 'Correlation Coefficients' 'value' 1 'tag' 'FCCA_CC' } {} ...
      { 'Style', 'text', 'string', 'General Setting', 'fontweight', 'bold' } {} {} ...                                       % 21
      { 'Style', 'text', 'string', '   Plot Status',  'FontAngle', 'italic' } ...                                            % 18
      { 'Style', 'checkbox', 'string' 'Plot All Smooth Curve' 'value' 0 'tag' 'AllSmoothPlot' } ...
      {'Style', 'checkbox', 'string' 'Cross Covariance and Cross Correlation' 'value' 0 'tag' 'S_A_H'} ... 
      };

      %%%% Calling (GUI)
   [ tmp1 tmp2 strhalt structout ] = inputgui( geometry, uilist, ...
           'pophelp(''pop_newtimef'');', 'ERP Analysis | Functional Canonical Correlation for EEGLAB - Version 0.1');

end  

      %%% Import GUI items into the variables.
      %%%%% General Settings

            nChannel = EEG.nbchan;           % Number of Channenls
            nTrials = EEG.trials;            % Number of Trials
            ChanLocs = EEG.chanlocs;         % Channels Location
            Events = EEG.event;              % Events
            Epochs = EEG.epoch;              % Epochs
            Times = EEG.times;               % EEG Times
            EData = EEG.data;                % EEG Dataset
            ICAWeights = EEG.icawinv;        % ICA weights
            ICAData = EEG.icaact;            % ICA Source activation

            %%% Input by User (General)
            %%% Channel (Electrode) Selection
            % K = 3;
            %%% Channel 1 
            SelectedChanIndex_1 = str2num(structout.SelectedChanIndex_1);
            SelectedChan_1 = ChanLocs(:,SelectedChanIndex_1) ;
            NSelectedChan_1 = length(SelectedChanIndex_1);
            
            %%% Channel 2
            SelectedChanIndex_2 = str2num(structout.SelectedChanIndex_2);
            SelectedChan_2 = ChanLocs(:,SelectedChanIndex_2) ;
            NSelectedChan_2 = length(SelectedChanIndex_2);

            %%% Event Selection
            % SelectedEventIndex = M;
            % SelectedEvent = Events(:,SelectedEventIndex);

            %%% Time Interval Selection 
            %%%% Channel 1
            TIS_1 = str2num(structout.TIS_1) ;
            TimeIntervalSele_1 = TIS_1; %%% Time Interval Selection (defualt: -200, 800 ms)
            TimeLO_1 = min(TimeIntervalSele_1) * 1.0e+03 ;
            TimeUP_1 = max(TimeIntervalSele_1) * 1.0e+03 ;
            TimeIntervalSeleInput_1      = Times(Times(1,:) >= TimeLO_1 & Times(1,:) <= TimeUP_1);
            TimeIntervalSeleIndexInput_1 = Times(1,:) >= TimeLO_1 & Times(1,:) <= TimeUP_1;

            %%%% Channel 2
            TIS_2 = str2num(structout.TIS_2) ;
            TimeIntervalSele_2 = TIS_2; %%% Time Interval Selection (defualt: -200, 800 ms)
            TimeLO_2 = min(TimeIntervalSele_2) * 1.0e+03 ;
            TimeUP_2 = max(TimeIntervalSele_2) * 1.0e+03 ;
            TimeIntervalSeleInput_2      = Times(Times(1,:) >= TimeLO_2 & Times(1,:) <= TimeUP_2);
            TimeIntervalSeleIndexInput_2 = Times(1,:) >= TimeLO_2 & Times(1,:) <= TimeUP_2;
            
            %%% Analysis Data Set with selected time intervals
            ADS_EData_1 = EData(:,TimeIntervalSeleIndexInput_1(1,:),:);
            ADS_EData_2 = EData(:,TimeIntervalSeleIndexInput_2(1,:),:);

            TypeBasisFunc_1 = structout.TypeBasisFunc_1;
            TypeBasisFunc_2 = structout.TypeBasisFunc_2;
            
            %%%% FCCA
            FCCA_STATUS = structout.FCCA_STATUS;  % FCCA done (1) or not (0);
            FCCA_NCAN = str2num(structout.FCCA_NCAN);      % Number of FCCA
            FCCA_NCAN_P = str2num(structout.FCCA_NCAN_P);      % FCCA Plot done (1) or not (0);
            FCCA_Diagonal = structout.FCCA_Diagonal;
            FCCA_offDiagonal = structout.FCCA_offDiagonal;
            FCCA_Diagonal_Scores  = structout.FCCA_Diagonal_Scores;
            FCCA_offDiagonal_Scores=  structout.FCCA_offDiagonal_Scores;
            FCCA_CC = structout.FCCA_CC;
            S_A_H = structout.S_A_H;

            
            %%% Analysis 
            %%%% Preprocessing
            %%% B-Spline
            %% Input
            %% Channel 1
            NumBasis_1 = str2num(structout.NB_1); % Number of Splines
            norder_1   = str2num(structout.NORD_1); % Number of Order
            TickGCV_1  = structout.TickGCV_1; % If GCV is ticked (TRUE-Default), the parameters are optimised with GCV, O.W. it is FALSE.
            PlotGCV_1 = structout.PlotGCV_1; % Plot GCV against Lambda Status (Deafult=0), If PlotGCV=1, it plots
            rng_1 = [min(TimeIntervalSeleInput_1),max(TimeIntervalSeleInput_1)];  % Range of B-Spline
            ADS_EData_Bspline_1_out_1=cell(NSelectedChan_1,1,nTrials);

            %% Channel 2
            NumBasis_2 = str2num(structout.NB_2); % Number of Splines
            norder_2   = str2num(structout.NORD_2); % Number of Order
            TickGCV_2  = structout.TickGCV_2; % If GCV is ticked (TRUE-Default), the parameters are optimised with GCV, O.W. it is FALSE.
            PlotGCV_2 = structout.PlotGCV_2; % Plot GCV against Lambda Status (Deafult=0), If PlotGCV=1, it plots
            rng_2 = [min(TimeIntervalSeleInput_2),max(TimeIntervalSeleInput_1)];  % Range of B-Spline
            ADS_EData_Bspline_1_out_2=cell(NSelectedChan_2,1,nTrials);

                               
            AllSmoothPlot = structout.AllSmoothPlot; % Plot all smoothed curves (default = 0 ) , if AllSmoothPlot=1, plot all smoothed curves. 
            %% Anlaysis Code
            %% for nch = 1:NSelectedChan
            %% 
            
                SelectedChnInd_1 = SelectedChanIndex_1(NSelectedChan_1);
                SelectedChanLabels_1 = SelectedChan_1(NSelectedChan_1).labels;
                ADS_EData_Bspline_1 = ADS_EData_1(SelectedChnInd_1,:,:);

                if TypeBasisFunc_1 == 1
                 %%% Defining B-Spline Basis  
                 %%% Preparing Data
                 wbasis_1 = create_bspline_basis(rng_1, NumBasis_1, norder_1);
                 ADS_EData_Bspline_1 = double(squeeze(ADS_EData_Bspline_1));
                 ADS_EData_Bspline_1 = ADS_EData_Bspline_1';
                 ADS_EData_Bspline_1 = double(ADS_EData_Bspline_1');
                 TimeIntervalSeleInput_BSpline_1 = double(TimeIntervalSeleInput_1');
                end 
                if TypeBasisFunc_1 == 2
                 %%% Defining Fourier Basis  
                 %%% Preparing Data
                wbasis_1 = create_fourier_basis(rng_1, NumBasis_1);
                 ADS_EData_Bspline_1 = double(squeeze(ADS_EData_Bspline_1));
                 ADS_EData_Bspline_1 = ADS_EData_Bspline_1';
                 ADS_EData_Bspline_1 = double(ADS_EData_Bspline_1');
                 TimeIntervalSeleInput_BSpline_1 = double(TimeIntervalSeleInput_1');
                end
               
                 if TickGCV_1 == 1
                     fprintf('The GCV is selected for the first channel.\n')

                           n_1 = length(ADS_EData_Bspline_1(1,:));
                      loglam_1 = (-9:0.25:9)'; %  set up the range of log lambda values
                        nlam_1 = length(loglam_1);
                         Lfd_1 = 4;

                    dfsave_1  = zeros(nlam_1,1);
                    gcvsave_1 = zeros(nlam_1,n_1)';
                    MSEsave_1 = zeros(nlam_1,n_1)';

                    %  loop through the log lambda values (Modify Later this part for parameters)
                    for ilam_1=1:nlam_1
                        lambda_1 = 10^loglam_1(ilam_1);
                        hgtfdPar_1 = fdPar(wbasis_1, Lfd_1, lambda_1);
                        [hgtfd_1, df_1, gcv_1, coef_1,SSE_1] = smooth_basis(TimeIntervalSeleInput_BSpline_1,ADS_EData_Bspline_1,hgtfdPar_1);
                        accest_1 = eval_fd(TimeIntervalSeleInput_BSpline_1, hgtfd_1, 2);
                        dfsave_1(ilam_1)    = df_1;
                        gcvsave_1(:,ilam_1)   = gcv_1;  % check later for mean 
                        MSEsave_1(:,ilam_1)   = SSE_1;  % check later for mean 
                        fprintf(['Channel Number 1 : ', num2str(SelectedChnInd_1),' , ', num2str(ilam_1),'-Log Lambda : ',num2str(loglam_1(ilam_1)), ' , Mean GCV of Epochs : ',num2str(round(mean(gcvsave_1(:,ilam_1)),4)),' , Mean SSE of Epochs: ',num2str(round(mean(MSEsave_1(:,ilam_1)),4))  ,'\n']);                       

                    end

                    [MVALVE_1,MININDEX_1] =  min(mean(gcvsave_1));
                    LOGLAMMIN_1 = loglam_1(MININDEX_1);

                    lambdaSelected_1  = 10^loglam_1(MININDEX_1);
                    hgtfdPar_1 = fdPar(wbasis_1, Lfd_1, lambdaSelected_1);
                    ADS_EData_Bspline_1_1 = smooth_basis(TimeIntervalSeleInput_BSpline_1,ADS_EData_Bspline_1,hgtfdPar_1);
                    fprintf(['*** Channel Number 1: ', num2str(SelectedChnInd_1),' , Minimum log Lambda 1: ',num2str(LOGLAMMIN_1), ' , Minimum GCV 1: ', num2str(MVALVE_1)]);

                    if PlotGCV_1 == 1
                    fprintf('The plot GCV is selected for the first channel.\n')
                    %%% Plot GCV (Modify)
                    figure(NSelectedChan_1)
                    plot(loglam_1,mean(gcvsave_1)' )        
                    axis([ (min(loglam_1)-3) (max(loglam_1)+3) (min(mean(gcvsave_1)-5)) (max(mean(gcvsave_1)+5)) ])
                    suptitle(['Channel Name: ',SelectedChanLabels_1])
                    title(['The minimum GCV (', num2str(MVALVE_1) ,') at Lambda:' num2str(loglam_1(MININDEX_1)) ])
                    yy_1 = MVALVE_1;
                    xx_1 = loglam_1(MININDEX_1);
                   line([xx_1, xx_1], [min(mean(gcvsave_1))-20 ,max(mean(gcvsave_1))+20],'Color','red','LineStyle','--')
                   line([min(loglam_1)-10, max(loglam_1)+10], [yy_1, yy_1],'Color','red','LineStyle','--')
                   ylabel('GCV') 
                   xlabel('Logarithm Lambda') 
                   end

                    %%% Plot MSE
                    %%%% plot(loglam,mean(MSEsave)) --> Based on Values
                    %%%% plot(mean(MSEsave))  %% --> Based on Index

                 elseif TickGCV_1 == 0
                      fprintf('The GCV is not selected for the first channel.\n')
                    %%% Smoothing Data
                    [ADS_EData_Bspline_1_1,df_1_1, gcv_1_1,coef_1_1,SSE_1_1] = smooth_basis(TimeIntervalSeleInput_BSpline_1,ADS_EData_Bspline_1,wbasis_1);
                 LOGLAMMIN_1 = 0;
                 end

                    %%% Output
                    EEGDatafd_names_1{1} = 'Time (ms)';
                    EEGDatafd_names_1{2} = SelectedChanLabels_1 ;
                    EEGDatafd_names_1{3} = '\mu. volt';
                    ADS_EData_Bspline_1_1 = putnames(ADS_EData_Bspline_1_1, EEGDatafd_names_1);

                                      %%% Output
                   %% Save all fd Object 
                   ADS_EData_Bspline_1_out_1(NSelectedChan_1,1,:) = fd2cell(ADS_EData_Bspline_1_1);


            %% Channel 2
           
                SelectedChnInd_2 = SelectedChanIndex_2(NSelectedChan_2);
                SelectedChanLabels_2 = SelectedChan_2(NSelectedChan_2).labels;
                ADS_EData_Bspline_2 = ADS_EData_2(SelectedChnInd_2,:,:);

                if TypeBasisFunc_2 == 1
                 %%% Defining B-Spline Basis  
                 %%% Preparing Data
                 wbasis_2 = create_bspline_basis(rng_2, NumBasis_2, norder_2);
                 ADS_EData_Bspline_2 = double(squeeze(ADS_EData_Bspline_2));
                 ADS_EData_Bspline_2 = ADS_EData_Bspline_2';
                 ADS_EData_Bspline_2 = double(ADS_EData_Bspline_2');
                 TimeIntervalSeleInput_BSpline_2 = double(TimeIntervalSeleInput_2'); 
                end 
                if TypeBasisFunc_2 == 2
                 %%% Defining Fourier Basis  
                 %%% Preparing Data
                 wbasis_2 = create_fourier_basis(rng_2, NumBasis_2);
                 ADS_EData_Bspline_2 = double(squeeze(ADS_EData_Bspline_2));
                 ADS_EData_Bspline_2 = ADS_EData_Bspline_2';
                 ADS_EData_Bspline_2 = double(ADS_EData_Bspline_2');
                 TimeIntervalSeleInput_BSpline_2 = double(TimeIntervalSeleInput_2'); 
                end


                 if TickGCV_2 == 1
                      fprintf('The GCV is selected for the second channel.\n')

                           n_2 = length(ADS_EData_Bspline_2(1,:));
                      loglam_2 = (-9:0.25:9)'; %  set up the range of log lambda values
                        nlam_2 = length(loglam_2);
                         Lfd_2 = 4;

                    dfsave_2  = zeros(nlam_2,1);
                    gcvsave_2 = zeros(nlam_2,n_2)';
                    MSEsave_2 = zeros(nlam_2,n_2)';

                    %  loop through the log lambda values (Modify Later this part for parameters)
                    for ilam_2=1:nlam_2
                        lambda_2 = 10^loglam_2(ilam_2);
                        hgtfdPar_2 = fdPar(wbasis_2, Lfd_2, lambda_2);
                        [hgtfd_2, df_2, gcv_2, coef_2,SSE_2] = smooth_basis(TimeIntervalSeleInput_BSpline_2,ADS_EData_Bspline_2,hgtfdPar_2);
                        accest_2 = eval_fd(TimeIntervalSeleInput_BSpline_2, hgtfd_2, 2);
                        dfsave_2(ilam_2)    = df_2;
                        gcvsave_2(:,ilam_2)   = gcv_2;  % check later for mean 
                        MSEsave_2(:,ilam_2)   = SSE_2;  % check later for mean 
                        fprintf(['Channel Number 2 : ', num2str(SelectedChnInd_2),' , ', num2str(ilam_2),'-Log Lambda : ',num2str(loglam_2(ilam_2)), ' , Mean GCV of Epochs : ',num2str(round(mean(gcvsave_2(:,ilam_2)),4)),' , Mean SSE of Epochs: ',num2str(round(mean(MSEsave_2(:,ilam_2)),4))  ,'\n']);                       
                    end

                    [MVALVE_2,MININDEX_2] =  min(mean(gcvsave_2));
                    LOGLAMMIN_2 = loglam_2(MININDEX_2);

                    lambdaSelected_2  = 10^loglam_2(MININDEX_2);
                    hgtfdPar_2 = fdPar(wbasis_2, Lfd_2, lambdaSelected_2);
                    ADS_EData_Bspline_1_2 = smooth_basis(TimeIntervalSeleInput_BSpline_2,ADS_EData_Bspline_2,hgtfdPar_2);
                    fprintf(['*** Channel Number 2: ', num2str(SelectedChnInd_2),' , Minimum log Lambda 2: ',num2str(LOGLAMMIN_2), ' , Minimum GCV 2: ', num2str(MVALVE_2)]);

                    
                    
                    if PlotGCV_2 == 1
                    %%% Plot GCV (Modify)
                    fprintf('The GCV plot is selected for the second channel.\n')
                    figure(NSelectedChan_2)
                    plot(loglam_2,mean(gcvsave_2)' )        
                    axis([ (min(loglam_2)-3) (max(loglam_2)+3) (min(mean(gcvsave_2)-5)) (max(mean(gcvsave_2)+5)) ])
                    suptitle(['Channel Name: ',SelectedChanLabels_2])
                    title(['The minimum GCV (', num2str(MVALVE_2) ,') at Lambda:' num2str(loglam_2(MININDEX_2)) ])
                    yy_2 = MVALVE_2;
                    xx_2 = loglam_2(MININDEX_2);
                   line([xx_2, xx_2], [min(mean(gcvsave_2))-20 ,max(mean(gcvsave_2))+20],'Color','red','LineStyle','--')
                   line([min(loglam_2)-10, max(loglam_2)+10], [yy_2, yy_2],'Color','red','LineStyle','--')
                   ylabel('GCV') 
                   xlabel('Logarithm Lambda') 
                   end

                    %%% Plot MSE
                    %%%% plot(loglam,mean(MSEsave)) --> Based on Values
                    %%%% plot(mean(MSEsave))  %% --> Based on Index

                 elseif TickGCV_2 == 0
                        fprintf('The GCV is not selected for second channel.\n')
                    %%% Smoothing Data
                    [ADS_EData_Bspline_1_2,df_1_2, gcv_1_2,coef_1_2,SSE_1_2] = smooth_basis(TimeIntervalSeleInput_BSpline_2,ADS_EData_Bspline_2,wbasis_2);
                 LOGLAMMIN_2 = 0;
                 end

                    %%% Output
                    EEGDatafd_names_2{1} = 'Time (ms)';
                    EEGDatafd_names_2{2} = SelectedChanLabels_2 ;
                    EEGDatafd_names_2{3} = '\mu. volt';
                    ADS_EData_Bspline_1_2 = putnames(ADS_EData_Bspline_1_2, EEGDatafd_names_2);

                   %% Save all fd Object 
                   ADS_EData_Bspline_1_out_2(NSelectedChan_2,1,:) = fd2cell(ADS_EData_Bspline_1_2);

                   
                   %%%% Plot each Curve Separatley  
                   %%% plotfit_fd(ADS_EData_Bspline, TimeIntervalSeleInput_BSpline, ADS_EData_Bspline_1)
                   if AllSmoothPlot == 1
                   fprintf('Plot all smoothed function is selected. \n');
                   %%% Curve 1     
                   AllSmoothPlot_1 = getcoef(ADS_EData_Bspline_1_1);
                   AllSmoothPlot_range_1 = getbasisrange(getbasis(ADS_EData_Bspline_1_1));
                   AllSmoothPlot_nbasis_1 = getnbasis(getbasis(ADS_EData_Bspline_1_1));
                   AllSmoothPlot_nx_1 = max([501, 10*AllSmoothPlot_nbasis_1+1]);
                   AllSmoothPlot_x_1  = linspace(AllSmoothPlot_range_1(1),AllSmoothPlot_range_1(2),AllSmoothPlot_nx_1)';
                   AllSmoothPlot_Lfdobj_1 = int2Lfd(int2Lfd(0));                           
                   AllSmoothPlot_fdmat_1 = eval_fd(AllSmoothPlot_x_1,ADS_EData_Bspline_1_1, AllSmoothPlot_Lfdobj_1);
                   AllSmoothPlot_max_1 = max(AllSmoothPlot_fdmat_1);
                   AllSmoothPlot_min_1 = min(AllSmoothPlot_fdmat_1);
                   
                   MaxS_1 = max(AllSmoothPlot_max_1);
                   MinS_1 = min(AllSmoothPlot_min_1);
                   
                   if TypeBasisFunc_1 == 1 
                    MName_1 = 'B-Spline';
                   end
                   
                   if TypeBasisFunc_1 == 2 
                    MName_1 = 'Fourier';
                   end

                   MinT_1 = min(TimeIntervalSeleInput_BSpline_1);
                   MaxT_1 = max(TimeIntervalSeleInput_BSpline_1);
                   
                  %%% Curve 2  
                   AllSmoothPlot_2 = getcoef(ADS_EData_Bspline_1_2);
                   AllSmoothPlot_range_2 = getbasisrange(getbasis(ADS_EData_Bspline_1_2));
                   AllSmoothPlot_nbasis_2 = getnbasis(getbasis(ADS_EData_Bspline_1_2));
                   AllSmoothPlot_nx_2 = max([501, 10*AllSmoothPlot_nbasis_2+1]);
                   AllSmoothPlot_x_2  = linspace(AllSmoothPlot_range_2(1),AllSmoothPlot_range_2(2),AllSmoothPlot_nx_2)';
                   AllSmoothPlot_Lfdobj_2 = int2Lfd(int2Lfd(0));                           
                   AllSmoothPlot_fdmat_2 = eval_fd(AllSmoothPlot_x_2,ADS_EData_Bspline_1_2, AllSmoothPlot_Lfdobj_2);
                   AllSmoothPlot_max_2 = max(AllSmoothPlot_fdmat_2);
                   AllSmoothPlot_min_2 = min(AllSmoothPlot_fdmat_2);
                   
                   MaxS_2 = max(AllSmoothPlot_max_2);
                   MinS_2 = min(AllSmoothPlot_min_2);
                   
                   if TypeBasisFunc_2 == 1 
                    MName_2 = 'B-Spline';
                   end
                   
                   if TypeBasisFunc_2 == 2 
                    MName_2 = 'Fourier';
                   end

                   MinT_2 = min(TimeIntervalSeleInput_BSpline_2);
                   MaxT_2 = max(TimeIntervalSeleInput_BSpline_2);
                   
                   figure;
                   plot(AllSmoothPlot_x_1 ,AllSmoothPlot_fdmat_1)
                   title(['The Smoothing Functions of ', SelectedChanLabels_1 ,' is:  ', MName_1 ] )
                   line([MinT_1,MaxT_1], [0,0],'Color',[0 0 0]+0.05,'LineStyle','--')
                   suptitle(['Channel Name: ',SelectedChanLabels_1,' (Best Log Lambda:' num2str(LOGLAMMIN_1),')'])
                   axis([MinT_1-10 MaxT_1+10 MinS_1-25 MaxS_1+25])
                   xlabel('Time (ms)') 
                   ylabel('Mean')      
                   
                   figure;
                   plot(AllSmoothPlot_x_2 ,AllSmoothPlot_fdmat_2)
                   title(['The Smoothing Functions of ', SelectedChanLabels_2 ,' is:  ', MName_2 ] )
                   line([MinT_2,MaxT_2], [0,0],'Color',[0 0 0]+0.05,'LineStyle','--')
                   suptitle(['Channel Name: ',SelectedChanLabels_2,' (Best Log Lambda:' num2str(LOGLAMMIN_2),')'])
                   axis([MinT_2-10 MaxT_2+10 MinS_2-25 MaxS_2+25])
                   xlabel('Time (ms)') 
                   ylabel('Mean') 
                   end

                   
                   
                   %%% FCCA
                   if FCCA_STATUS == 1
                   fprintf('The functional canonincal correlation is selected. \n');
                   FCCA_NCAN = FCCA_NCAN;
                   ccafdPar_1 = fdPar(wbasis_1, 2, 5e6);
                   ccafdPar_2 = fdPar(wbasis_2, 2, 5e6);

                   FCCA_IBG = cca_fd(ADS_EData_Bspline_1_1,ADS_EData_Bspline_1_2,FCCA_NCAN,ccafdPar_1,ccafdPar_2);

                  if FCCA_Diagonal == 1 
                  fprintf('Plot diagonal weight functions. \n');
                  for ii = 1:FCCA_NCAN_P
                       for jj = 1:FCCA_NCAN_P
                           if ii == jj         
                               fprintf(['Canonical weight function of ', num2str(ii),' against ',num2str(jj) , '.\n'])  
                               curve_ii = getcoef(FCCA_IBG.wtfdx(ii));
                               curve_jj = getcoef(FCCA_IBG.wtfdy(jj));
                               curve_ii_range = getbasisrange(getbasis(FCCA_IBG.wtfdx(ii)));
                               curve_jj_range = getbasisrange(getbasis(FCCA_IBG.wtfdy(jj)));
                               curve_ii_nbasis   = getnbasis(getbasis(FCCA_IBG.wtfdx(ii)));
                               curve_jj_nbasis   = getnbasis(getbasis(FCCA_IBG.wtfdy(jj)));
                               curve_ii_nx = max([501, 10*curve_ii_nbasis+1]);
                               curve_jj_nx = max([501, 10*curve_jj_nbasis+1]);
                               curve_ii_x     = linspace(curve_ii_range(1),curve_ii_range(2),curve_ii_nx)';
                               curve_ii_Lfdobj = int2Lfd(int2Lfd(0));                           
                               curve_ii_fdmat = eval_fd(curve_ii_x, FCCA_IBG.wtfdx(ii), curve_ii_Lfdobj);
                               curve_jj_x     = linspace(curve_jj_range(1),curve_jj_range(2),curve_jj_nx)';
                               curve_jj_Lfdobj = int2Lfd(int2Lfd(0));                           
                               curve_jj_fdmat = eval_fd(curve_jj_x, FCCA_IBG.wtfdy(jj), curve_jj_Lfdobj);
                                max_ii = max(curve_ii);
                                max_jj = max(curve_jj);
                                min_ii = min(curve_ii);
                                min_jj = min(curve_jj);
                                MinS = min(min_ii, min_jj);
                                MaxS = max(max_ii, max_jj);
                                MinT = min(min(TimeIntervalSeleInput_1),min(TimeIntervalSeleInput_2));
                                MaxT = max(max(TimeIntervalSeleInput_1),max(TimeIntervalSeleInput_2));
                                figure()
                                plot(curve_ii_x,curve_ii_fdmat,'color','black','LineWidth',1.5,'LineStyle' ,'-')
                                hold on 
                                plot(curve_jj_x,curve_jj_fdmat,'color','red','LineWidth',1.5,'LineStyle' ,'-')
                                title(['The ', num2str(ii) ,' Canonical weight function of ', SelectedChanLabels_1 ,' against the ', num2str(jj),' of ',SelectedChanLabels_2])
                                line([MinT,MaxT], [0,0],'Color',[0 0 0]+0.05,'LineStyle','--')
                                axis([MinT-10 MaxT+10 MinS-0.10 MaxS+0.10])
                                legend({[num2str(ii) ,' CWF ', SelectedChanLabels_1],[num2str(jj) ,' CWF ' , SelectedChanLabels_2]})
                                xlabel('Time (ms)') 
                                ylabel('Canonical Weight Function')
                               
                                                      
                           end 
                       end
                  end
              end 
                  

               if FCCA_offDiagonal == 1 
                  fprintf('Plot all weight functions. \n');
                  for ii = 1:FCCA_NCAN_P
                       for jj = 1:FCCA_NCAN_P
                          % if ii < jj  
                               fprintf(['Canonical weight function of ', num2str(ii),' against ',num2str(jj) , '.\n'])  
                               curve_ii = getcoef(FCCA_IBG.wtfdx(ii));
                               curve_jj = getcoef(FCCA_IBG.wtfdy(jj));

                               curve_ii_range = getbasisrange(getbasis(FCCA_IBG.wtfdx(ii)));
                               curve_jj_range = getbasisrange(getbasis(FCCA_IBG.wtfdy(jj)));

                               curve_ii_nbasis   = getnbasis(getbasis(FCCA_IBG.wtfdx(ii)));
                               curve_jj_nbasis   = getnbasis(getbasis(FCCA_IBG.wtfdy(jj)));

                               curve_ii_nx = max([501, 10*curve_ii_nbasis+1]);
                               curve_jj_nx = max([501, 10*curve_jj_nbasis+1]);

                               curve_ii_x     = linspace(curve_ii_range(1),curve_ii_range(2),curve_ii_nx)';
                               curve_ii_Lfdobj = int2Lfd(int2Lfd(0));                           
                               curve_ii_fdmat = eval_fd(curve_ii_x, FCCA_IBG.wtfdx(ii), curve_ii_Lfdobj);

                               curve_jj_x     = linspace(curve_jj_range(1),curve_jj_range(2),curve_jj_nx)';
                               curve_jj_Lfdobj = int2Lfd(int2Lfd(0));                           
                               curve_jj_fdmat = eval_fd(curve_jj_x, FCCA_IBG.wtfdy(jj), curve_jj_Lfdobj);
                           
                                max_ii = max(curve_ii);
                                max_jj = max(curve_jj);
                                min_ii = min(curve_ii);
                                min_jj = min(curve_jj);
                                
                                MinS = min(min_ii, min_jj);
                                MaxS = max(max_ii, max_jj);
                                MinT = min(min(TimeIntervalSeleInput_1),min(TimeIntervalSeleInput_2));
                                MaxT = max(max(TimeIntervalSeleInput_1),max(TimeIntervalSeleInput_2));
                                
                                figure()
                                plot(curve_ii_x,curve_ii_fdmat,'color','black','LineWidth',1.5,'LineStyle' ,'-')
                                hold on 
                                plot(curve_jj_x,curve_jj_fdmat,'color','red','LineWidth',1.5,'LineStyle' ,'-')

                                title(['The ', num2str(ii) ,' Canonical weight function of ', SelectedChanLabels_1 ,' against the ', num2str(jj),' of ',SelectedChanLabels_2])
                                line([MinT,MaxT], [0,0],'Color',[0 0 0]+0.05,'LineStyle','--')
                                axis([MinT-10 MaxT+10 MinS-0.10 MaxS+0.10])
                                legend({[num2str(ii) ,' CWF ', SelectedChanLabels_1],[num2str(jj) ,' CWF ' , SelectedChanLabels_2]})
                                xlabel('Time (ms)') 
                                ylabel('Canonical Weight Function')

                                
                                % end 
                       end
                  end
               end 
               
               
               if FCCA_Diagonal_Scores == 1 
                  fprintf('Diagonal Canonical scores function is selected .\n')  
                  for ii = 1:FCCA_NCAN_P
                       for jj = 1:FCCA_NCAN_P
                           if ii == jj   
                                strName = num2str((1:size(FCCA_IBG.varx(:,ii)))');
                                fprintf(['Canonical scores function ', num2str(ii),' against ',num2str(jj) , '.\n'])  
                                figure()
                                plot(FCCA_IBG.varx(:,ii),FCCA_IBG.vary(:,jj),".","MarkerSize",14,'color','black')
                                line([0, 0], [min(FCCA_IBG.vary(:,jj))-300,max(FCCA_IBG.vary(:,jj))+300],'Color','red','LineStyle','--')
                                line([min(FCCA_IBG.varx(:,ii))-300,max(FCCA_IBG.varx(:,ii))+300], [0,0],'Color','red','LineStyle','--')
                                title(['The ', num2str(ii) ,' Canonical scores function of ', SelectedChanLabels_1 ,' against the ', num2str(jj),' of ',SelectedChanLabels_2])
                                xlabel(['Scores ', num2str(ii), ' of ',SelectedChanLabels_1]) 
                                ylabel(['Scores ', num2str(jj), ' of ',SelectedChanLabels_2])
                                axis([ (min(FCCA_IBG.varx(:,ii))-100) (max(FCCA_IBG.varx(:,ii))+100) (min(FCCA_IBG.vary(:,jj))-100) (max(FCCA_IBG.vary(:,jj))+100) ])
                                text(FCCA_IBG.varx(:,ii),FCCA_IBG.vary(:,jj),strName,'Color', 'blue','FontSize',8) 

                           end 
                       end
                  end
              end 
                  

               if FCCA_offDiagonal_Scores == 1 
                  fprintf('All Canonical scores function is selected .\n')  
                  for ii = 1:FCCA_NCAN_P
                       for jj = 1:FCCA_NCAN_P
                         %  if ii < jj
                                fprintf(['Canonical scores function ', num2str(ii),' against ',num2str(jj) , '.\n'])  
                                strName = num2str((1:size(FCCA_IBG.varx(:,ii)))');
                                figure()
                                plot(FCCA_IBG.varx(:,ii),FCCA_IBG.vary(:,jj),".","MarkerSize",14,'color','black')
                                line([0, 0], [min(FCCA_IBG.vary(:,jj))-300,max(FCCA_IBG.vary(:,jj))+300],'Color','red','LineStyle','--')
                                line([min(FCCA_IBG.varx(:,ii))-300,max(FCCA_IBG.varx(:,ii))+300], [0,0],'Color','red','LineStyle','--')
                                title(['The ', num2str(ii) ,' Canonical scores function of ', SelectedChanLabels_1 ,' against the ', num2str(jj),' of ',SelectedChanLabels_2])
                                xlabel(['Scores ', num2str(ii), ' of ',SelectedChanLabels_1]) 
                                ylabel(['Scores ', num2str(jj), ' of ',SelectedChanLabels_2])
                                axis([ (min(FCCA_IBG.varx(:,ii))-100) (max(FCCA_IBG.varx(:,ii))+100) (min(FCCA_IBG.vary(:,jj))-100) (max(FCCA_IBG.vary(:,jj))+100) ])
                                text(FCCA_IBG.varx(:,ii),FCCA_IBG.vary(:,jj),strName,'Color', 'blue','FontSize',8) 
                            % end 
                       end
                  end
               end 
               
               
               
               if FCCA_CC == 1 
                   fprintf('The correlation coefficients plot is selected.\n')  
                   SZM = size(FCCA_IBG.corrs);
                   rpt =   [(1:SZM)',round(FCCA_IBG.corrs*100,2)];
                   tableReport = array2table(rpt,'VariableNames',{'Basis_Functions','Canonical_Correlation'})
                             
                   figure()
                   plot(FCCA_IBG.corrs,"+-",'Color','black','LineWidth',2)
                   line([0,SZM(1)], [0.9 0.9],'Color','blue','LineStyle','--','LineWidth',2)
                   line([sum(FCCA_IBG.corrs >= 0.9),sum(FCCA_IBG.corrs >= 0.9)], [0 1],'Color','blue','LineStyle','--','LineWidth',2)
                  % line([0,SZM(1)], [0.7 0.7],'Color','green','LineStyle','--','LineWidth',2)
                  % line([sum(FCCA_IBG.corrs >= 0.7),sum(FCCA_IBG.corrs >= 0.8)], [0 1],'Color','green','LineStyle','--','LineWidth',2)
                  % line([0,SZM(1)], [0.5 0.5],'Color','red','LineStyle','--','LineWidth',2)
                  % line([sum(FCCA_IBG.corrs >= 0.5),sum(FCCA_IBG.corrs >= 0.5)], [0 1],'Color','red','LineStyle','--','LineWidth',2)
                   title('The Canonical Correlation Coefficients')
                   xlabel('Number of Basis Function') 
                   ylabel('Correlation Coefficient') 
                   grid on
                   grid minor
               end 
      end 
         
                   
                if S_A_H==1 
                   fprintf('The cross-covariance and cross-correlation is selected.\n')  
                   fprintf(['Channel ' ,SelectedChanLabels_1 , ' vs channel ', SelectedChanLabels_2 ,'.\n']) 
                   curve_ii = getcoef(ADS_EData_Bspline_1_1);
                   curve_jj = getcoef(ADS_EData_Bspline_1_2);
                   curve_ii_range = getbasisrange(getbasis(ADS_EData_Bspline_1_1));
                   curve_jj_range = getbasisrange(getbasis(ADS_EData_Bspline_1_2));
                   curve_ii_nbasis   = getnbasis(getbasis(ADS_EData_Bspline_1_1));
                   curve_jj_nbasis   = getnbasis(getbasis(ADS_EData_Bspline_1_2));
                   curve_ii_nx = max([501, 10*curve_ii_nbasis+1]);
                   curve_jj_nx = max([501, 10*curve_jj_nbasis+1]);
                   curve_ii_x     = linspace(curve_ii_range(1),curve_ii_range(2),curve_ii_nx)';
                   curve_ii_Lfdobj = int2Lfd(int2Lfd(0));                           
                   curve_jj_x     = linspace(curve_jj_range(1),curve_jj_range(2),curve_jj_nx)';
                   %%% Cross-Covariance (cor_fd)
                   CC_1 = var_fd(ADS_EData_Bspline_1_1,ADS_EData_Bspline_1_2);
                   CC_Result = eval_bifd(curve_ii_x,curve_jj_x, CC_1);
                   CC_2 = cor_fd(curve_ii_x,ADS_EData_Bspline_1_1,curve_jj_x,ADS_EData_Bspline_1_2);
                   figure();
                   subplot(2,2,1);
                   surf(curve_ii_x, curve_jj_x, CC_Result,'EdgeColor','interp','FaceColor','flat','FaceLighting','gouraud','FaceAlpha',0.5);
                   colorbar;
                   caxis([min(min(CC_Result)) max(max(CC_Result))])
                   title(['The Cross-Covariance Surface between ',SelectedChanLabels_1,' and ',SelectedChanLabels_2,'.'])
                   subplot(2,2,2);
                   contour(curve_ii_x, curve_jj_x, CC_Result);
                   colorbar
                   caxis([min(min(CC_Result)) max(max(CC_Result))])
                   title(['The Cross-Covariance Image between ',SelectedChanLabels_1,' and ',SelectedChanLabels_2,'.'])
                   subplot(2,2,3);
                   surf(CC_2,'EdgeColor','interp','FaceColor','flat','FaceLighting','gouraud','FaceAlpha',0.5);
                   colorbar;
                   caxis([-1 1])
                   title(['The Cross-Correlation Surface between ',SelectedChanLabels_1,' and ',SelectedChanLabels_2,'.'])
                   subplot(2,2,4);
                   contour(CC_2);
                   colorbar
                   caxis([-1 1])
                   title(['The Cross-Correlation Image between ',SelectedChanLabels_1,' and ',SelectedChanLabels_2,'.'])

               end 
end 