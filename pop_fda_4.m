% pop_fda_4 - perform ERP analysis with functional data analysis (FDA).  
%             It also smooth curves with B-Spline and Fourier basis functions. 
%             The estimation of penalty term is obtained with generalized cross
%             validavtion (GCV). The derivatives and phase-plane plot are estimated. 
%             It also computes the functional principal component analysis (FPCA) on 
%             channels. 
%
% Usage:
%  pop_fda_4(EEG); % pop up window asking users to select method
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
%   pop_fda_4(EEG); % pop up window asking users to select method
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

function com = pop_fda_4(EEG, varargin)
%%%% Tools -> ERP Analysis -> Descriptive

fprintf('ERP Analysis \n');
fprintf('Descriptive \n');


com = '';
if nargin < 1
    help pop_sourcereconstruction;
    return
end

if nargin < 2
    %%%% GUI 
    g = [1 0.5 0.5];
    geometry = { g g g g g g g g g g g g g g g g g g g g g};
    uilist = { ...
      { 'Style', 'text', 'string', 'Channel Selection', 'fontweight', 'bold' } {} {} ...                                              %1
      { 'Style', 'text', 'string', '   Channel number (At least two channels)', 'FontAngle', 'italic'  } ...                          %2
      { 'Style', 'edit', 'string', '1:10' 'tag' 'SelectedChanIndex'} ...
      {'Style', 'checkbox', 'string' 'All Channels' 'value' 1 'tag' 'AllChannels'} ...
      ...
      { 'Style', 'text', 'string', '   Time limits [min max] (msec)', 'FontAngle', 'italic' } ...                                     %3
      { 'Style', 'edit', 'string', '-0.2 0.800 ' 'tag' 'TIS' } {} ...
      ...
      { 'Style', 'text', 'string', 'Preprocessing (Step 1) ', 'fontweight', 'bold' } {} {} ...                                        %4
      { 'Style', 'text', 'string', '   Choose Basis Function (Default: B-Spline)', 'FontAngle', 'italic' } ...                        %5
      { 'Style', 'popupmenu', 'string', 'B-Spline|Fourier' , 'tag' 'TypeBasisFunc' } {} ...
      { 'Style', 'text', 'string', '   Number of Basis', 'FontAngle', 'italic' } ...                                                  %6
      { 'Style', 'edit', 'string', '120' 'tag' 'NB'  } {} ...
      ...
      { 'Style', 'text', 'string', '   Order of Basis','FontAngle', 'italic' } ...                                                    %7
      { 'Style', 'edit', 'string', '6' 'tag' 'NORD' } {} ...
      { 'Style', 'text', 'string', '   Parameter Estimation', 'FontAngle', 'italic' } ...                                             %8
      { 'Style', 'checkbox', 'string' 'GCV' 'value' 1 'tag' 'TickGCV' } ...
      { 'Style', 'checkbox', 'string' 'Plot GCV' 'value' 1 'tag' 'PlotGCV' } ...
      { 'Style', 'text', 'string', 'Descriptive Analysis', 'fontweight', 'bold' } {} {} ...                                           %9
      { 'Style', 'text', 'string', '   ERP Plot', 'FontAngle', 'italic' } ...                                                         %10
      { 'Style', 'checkbox', 'string' 'Mean and Standard Deviation ' 'value' 0 'tag' 'MP' } {} ...
      { 'Style', 'text', 'string', '   Status', 'FontAngle', 'italic' } ...                                                           %11
      { 'Style', 'checkbox', 'string' 'Derivative and Phase-Plane Plot' 'value' 1 'tag' 'Deri_PP_plot' } {} ...
      { 'Style', 'text', 'string', '   Derivative (order)', 'FontAngle', 'italic' } ...                                               %12
      { 'Style', 'edit', 'string', '1' 'tag' 'NODERIVATIVE' } {} ...
      { 'Style', 'text', 'string', '   ERP Derivative', 'FontAngle', 'italic' } ...                                                   %13
      {'Style', 'checkbox', 'string' 'Plot' 'value' 0 'tag' 'PPL_DER'} {} ...
      { 'Style', 'text', 'string', '   ERP Phase-Plane', 'FontAngle', 'italic' } ...                                                  %14
      { 'Style', 'checkbox', 'string' 'Plot' 'value' 0 'tag' 'PPL' } {} ... 
      { 'Style', 'text', 'string', 'Functional Principal Components Analysis (FPCA)', 'fontweight', 'bold' } {} {} ...                %15
      { 'Style', 'text', 'string', '   FPCA Status',  'FontAngle', 'italic' } ...                                                     %16
      { 'Style', 'checkbox', 'string' 'Done' 'value' 1 'tag' 'FPCA_STATUS' } {} ...
      { 'Style', 'text', 'string', '   Number of FPCA', 'FontAngle', 'italic' } ...                                                   %17
      { 'Style', 'edit', 'string', '4' 'tag' 'FPCA_nFPCA'} {} ...
      { 'Style', 'text', 'string', '   FPCA Descriptive', 'FontAngle', 'italic' } ...                                                 %18
      { 'Style', 'checkbox', 'string' 'Plot Scree plot' 'value' 1 'tag' 'FPCA_PlotFPCA_SCREE' } ...
      { 'Style', 'checkbox', 'string' 'Plot Eigefunction' 'value' 0 'tag' 'FPCA_PlotFPCA' } ...
      { 'Style', 'text', 'string', '   Rotations', 'FontAngle', 'italic' } ...                                                        %19
      { 'Style', 'checkbox', 'string' 'Varimax Rotation' 'value' 1 'tag' 'FPCA_PlotFPCA_VARMAX' }  ...
      { 'Style', 'checkbox', 'string' 'Plot Scree plot' 'value' 0 'tag' 'FPCA_VARIMAX_PlotFPCA_SCREE' } ...      
      { 'Style', 'text', 'string', '   Plot Scores', 'FontAngle', 'italic' } ...                                                      %20
      { 'Style', 'checkbox', 'string' 'Scores' 'value' 0 'tag' 'FPCA_Scores' } ...
      { 'Style', 'checkbox', 'string' 'Varimax Scores' 'value' 0 'tag' 'FPCA_Scores_VARMAX'} ...   
      { 'Style', 'text', 'string', '   Plot 3D Scores (First Three FPCs)', 'FontAngle', 'italic' } ...                                                      %20
      { 'Style', 'checkbox', 'string' 'Scores' 'value' 1 'tag' 'FPCA_Scores_3D' } ...
      { 'Style', 'checkbox', 'string' 'Varimax Scores' 'value' 0 'tag' 'FPCA_Scores_VARMAX_3D'} ...  
      };

      %%%% Calling (GUI)
      [ tmp1 tmp2 strhalt structout ] = inputgui( geometry, uilist, ...
           'pophelp(''pop_newtimef'');', 'ERP Analysis | Descriptive of Functional Data Analysis for EEGLAB - Version 0.1');
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
            
      %%%%% User's entries 
      %%%%%%%% Channel (Electrode) Selection
            AllChannels = structout.AllChannels;       % Consider All Channels      (Default: 1)
            if AllChannels == 1 
                NSelectedChan = nChannel;
                SelectedChanIndex = 1:nChannel;
                SelectedChan = ChanLocs(:,SelectedChanIndex) ;
            end
            if  AllChannels == 0
                SelectedChanIndex = str2num(structout.SelectedChanIndex);
                SelectedChan = ChanLocs(:,SelectedChanIndex) ;
                NSelectedChan = length(SelectedChanIndex);
            end 
           
      %%%%%%%% Time Interval Selection
            TIS = str2num(structout.TIS) ;
            TimeIntervalSele = TIS;       %%% Time Interval Selection (defualt: -200, 800 ms)
            TimeLO = min(TimeIntervalSele) * 1.0e+03 ;
            TimeUP = max(TimeIntervalSele) * 1.0e+03 ;
            TimeIntervalSeleInput      = Times(Times(1,:) >= TimeLO & Times(1,:) <= TimeUP);
            TimeIntervalSeleIndexInput = Times(1,:) >= TimeLO & Times(1,:) <= TimeUP;
            ADS_EData = EData(:,TimeIntervalSeleIndexInput(1,:),:);   %%% Analysis Data Set with selected time intervals
      %%% Analysis
      %%%%%%%% Inputs Variables
            NumBasis = str2num(structout.NB);                               % Number of Basis Functions (Default: 120)
            norder   = str2num(structout.NORD);                             % Order of Basis Functions (Default: 6)
            TypeBasisFunc = structout.TypeBasisFunc;                        % Type of Basis Function     (Default: B-Spline)
            TickGCV  = structout.TickGCV;                                   % Generalised Cross Validation (GCV) status for estimation the penalty term (Lambda): (Default:1->Done)
            PlotGCV = structout.PlotGCV;                                    % Plot GCV against log Lambda Status (Default:1->Don)
            rng = [min(TimeIntervalSeleInput),max(TimeIntervalSeleInput)];  % Range for basis function 
            ADS_EData_Bspline_1_out=cell(NSelectedChan,1,nTrials);          % output defintion  
            PPL= structout.PPL;                                             % Plot the Phase-Plane Plots (Default: 0)
            MP = structout.MP;                                              % Plot Mean Function         (Default: 0)
            NODERIVATIVE = str2num(structout.NODERIVATIVE);                 % The order of Derivaive     (Default: 1)
            PPL_DER = structout.PPL_DER;                                    % Plot the Derivative        (Default: 0)
            Deri_PP_plot = structout.Deri_PP_plot;                          % Calculate the Derivative and Phase-Plane Plots  (Default: 0) 
            FPCA_STATUS = structout.FPCA_STATUS;                            % FPCA Calculations          (Default: 0)
            FPCA_nFPCA = str2num(structout.FPCA_nFPCA);                              % Number of FPCA             (Default: 0)
            FPCA_PlotFPCA_SCREE = structout.FPCA_PlotFPCA_SCREE;            % Screeplot for FPCA         (Default: 0)
            FPCA_PlotFPCA = structout.FPCA_PlotFPCA;                        % Plot eigenfunction         (Default: 0)
            FPCA_PlotFPCA_VARMAX = structout.FPCA_PlotFPCA_VARMAX;          % Plot eigenfunction Varimax FPCA   (Default: 0)
            FPCA_VARIMAX_PlotFPCA_SCREE = structout.FPCA_VARIMAX_PlotFPCA_SCREE;   % Screeplot for Varimax FPCA   (Default: 0)
            FPCA_Scores = structout.FPCA_Scores;                            % Plot Scores FPCA
            FPCA_Scores_VARMAX = structout.FPCA_Scores_VARMAX;              % Plot Scores FPCA Varimax
            FPCA_Scores_3D        = structout.FPCA_Scores_3D;               % Plot FPCA Scores 3D 
            FPCA_Scores_VARMAX_3D = structout.FPCA_Scores_VARMAX_3D;        % Plot FPCA Varimax Scores 3D 
            if TypeBasisFunc == 1
               NameBasisFunc = "B-Spline";
             end 
             if TypeBasisFunc == 2
                NameBasisFunc = "Fourier";
             end 
      %%%%%%%% Start Analysis
            for nch = 1:NSelectedChan
      %%%%%%%%%%%%%%%% Selecting Channels
                SelectedChnInd = SelectedChanIndex(nch);
                SelectedChanLabels = SelectedChan(nch).labels;
                ADS_EData_Bspline = ADS_EData(SelectedChnInd,:,:);
      %%%%%%%%%%%%%%%% Basis Definition 
                if TypeBasisFunc == 1
      %%%%%%%%%%%%%%%%%%% B-Spline Basis   
                 wbasis = create_bspline_basis(rng, NumBasis, norder);
                 ADS_EData_Bspline = double(squeeze(ADS_EData_Bspline));
                 ADS_EData_Bspline = ADS_EData_Bspline';
                 ADS_EData_Bspline = double(ADS_EData_Bspline');
                 TimeIntervalSeleInput_BSpline = double(TimeIntervalSeleInput');
                end 
                if TypeBasisFunc == 2
      %%%%%%%%%%%%%%%%%%% Fourier Basis   
                 wbasis = create_fourier_basis(rng, NumBasis);
                 ADS_EData_Bspline = double(squeeze(ADS_EData_Bspline));
                 ADS_EData_Bspline = ADS_EData_Bspline';
                 ADS_EData_Bspline = double(ADS_EData_Bspline');
                 TimeIntervalSeleInput_BSpline = double(TimeIntervalSeleInput');
                end
      %%%%%%%%%%%%%%%% GCV for Lambda estimation 
                 if TickGCV == 1
                           n = length(ADS_EData_Bspline(1,:));
                      loglam = (-9:0.25:9)';
                        nlam = length(loglam);
                         Lfd = 4;
                    dfsave  = zeros(nlam,1);
                    gcvsave = zeros(nlam,n)';
                    MSEsave = zeros(nlam,n)';
                    for ilam=1:nlam
                        lambda = 10^loglam(ilam);
                        hgtfdPar = fdPar(wbasis, Lfd, lambda);
                        [hgtfd, df, gcv, coef,SSE] = smooth_basis(TimeIntervalSeleInput_BSpline,ADS_EData_Bspline,hgtfdPar);
                        accest = eval_fd(TimeIntervalSeleInput_BSpline, hgtfd, 2);
                        dfsave(ilam)    = df;
                        gcvsave(:,ilam)   = gcv; 
                        MSEsave(:,ilam)   = SSE;  
                        fprintf(['Channel Number : ', num2str(SelectedChnInd),' , ', num2str(ilam),'-Log Lambda : ',num2str(loglam(ilam)), ' , Mean GCV of Epochs : ',num2str(round(mean(gcvsave(:,ilam)),4)),' , Mean SSE of Epochs: ',num2str(round(mean(MSEsave(:,ilam)),4))  ,'\n']);
                    end
                    [MVALVE,MININDEX] =  min(mean(gcvsave));
                    LOGLAMMIN = loglam(MININDEX);
                    lambdaSelected  = 10^loglam(MININDEX);         % Minimum Lambda
                    hgtfdPar = fdPar(wbasis, Lfd, lambdaSelected);
                    ADS_EData_Bspline_1 = smooth_basis(TimeIntervalSeleInput_BSpline,ADS_EData_Bspline,hgtfdPar);
                    fprintf(['*** Channel Number: ', num2str(SelectedChnInd),' , Minimum log Lambda: ',num2str(LOGLAMMIN), ' , Minimum GCV: ', num2str(MVALVE)]);
      %%%%%%%%%%%%%%%%%%% Plot GCV against Log-Lambda   
                    if PlotGCV == 1
                    fprintf('The GCV Plot is selected.\n')
                    figure()
                    plot(loglam,mean(gcvsave)' )        
                    axis([ (min(loglam)-3) (max(loglam)+3) (min(mean(gcvsave)-5)) (max(mean(gcvsave)+5)) ])
                    suptitle(['Channel Name: ',SelectedChanLabels])
                    title(['The minimum GCV (', num2str(MVALVE) ,') at Lambda:' num2str(loglam(MININDEX)) ])
                    yy = MVALVE;
                    xx = loglam(MININDEX);
                   line([xx, xx], [min(mean(gcvsave))-20 ,max(mean(gcvsave))+20],'Color','red','LineStyle','--')
                   line([min(loglam)-10, max(loglam)+10], [yy, yy],'Color','red','LineStyle','--')
                   ylabel('GCV') 
                   xlabel('Logarithm Lambda') 
                    end
      %%%%%%%%%%%%%%%% not used GCV for Lambda estimation 
                elseif TickGCV == 0
                     fprintf('The GCV Plot is not selected.')
                    [ADS_EData_Bspline_1,df_1, gcv_1,coef_1,SSE_1] = smooth_basis(TimeIntervalSeleInput_BSpline,ADS_EData_Bspline,wbasis);
                     LOGLAMMIN = 0;
                 end
      %%%%%%%%%%%%%%%% Defining smoothed data
                    EEGDatafd_names{1} = 'Time (ms)';
                    EEGDatafd_names{2} = SelectedChanLabels ;
                    EEGDatafd_names{3} = '\mu. volt';
                    ADS_EData_Bspline_1 = putnames(ADS_EData_Bspline_1, EEGDatafd_names);

      %%%%%%%%%%%%%%% Mean of all smoothed epoched of a channel
                   Mean_FD = mean(ADS_EData_Bspline_1);
                   curve_Mean_FD = getcoef(Mean_FD);
                   curve_Mean_FD_range = getbasisrange(getbasis(Mean_FD));
                   curve_Mean_FD_nbasis = getnbasis(getbasis(Mean_FD));
                   curve_Mean_FD_nx = max([501, 10*curve_Mean_FD_nbasis+1]);
                   curve_Mean_FD_x  = linspace(curve_Mean_FD_range(1),curve_Mean_FD_range(2),curve_Mean_FD_nx)';
                   curve_Mean_FD_Lfdobj = int2Lfd(int2Lfd(0));                           
                   curve_Mean_FD_fdmat = eval_fd(curve_Mean_FD_x,Mean_FD, curve_Mean_FD_Lfdobj);
                   
      %%%%%%%%%%%%%%% STD of all smoothed epoched of a channel

                   STD_FD = std_fd(ADS_EData_Bspline_1);
                   curve_STD_FD = getcoef(STD_FD);
                   curve_STD_FD_range = getbasisrange(getbasis(STD_FD));
                   curve_STD_FD_nbasis = getnbasis(getbasis(STD_FD));
                   curve_STD_FD_nx = max([501, 10*curve_STD_FD_nbasis+1]);
                   curve_STD_FD_x  = linspace(curve_STD_FD_range(1),curve_STD_FD_range(2),curve_STD_FD_nx)';
                   curve_STD_FD_Lfdobj = int2Lfd(int2Lfd(0));                           
                   curve_STD_FD_fdmat = eval_fd(curve_STD_FD_x,STD_FD, curve_STD_FD_Lfdobj);

                   if nch == 1
                         fprintf('The number of channel is 1. \n')
                         TimePoints = size(curve_Mean_FD_x);
                         Mean_FD_Matrix = zeros( TimePoints(1), NSelectedChan+1);
                         Mean_FD_Matrix(:,1)  = curve_Mean_FD_x;
                         Mean_FD_Matrix(:,nch+1)  = curve_Mean_FD_fdmat;
                         
                         TimePoints = size(curve_STD_FD_x);
                         STD_FD_Matrix = zeros( TimePoints(1), NSelectedChan+1);
                         STD_FD_Matrix(:,1)  = curve_STD_FD_x;
                         STD_FD_Matrix(:,nch+1)  = curve_STD_FD_fdmat;
                   end 
                   
                   if nch > 1
                         fprintf('The number of channel is greater than 1. \n')
                         Mean_FD_Matrix(:,nch+1)  = curve_Mean_FD_fdmat;
                         STD_FD_Matrix(:,nch+1)   = curve_STD_FD_fdmat;
                   end
            end 
      %%%%%%%%%%%%%%%% Calculate and Plot mean and standard deviation function of smoothed curve 
                   if MP == 1
                      if NSelectedChan == 1
                        fprintf('The number of channel is 1.\n')  
                        fprintf('The functioal mean and standard deviation are plotted.\n')  
                        figure()
                        subplot(2,2,1)
                        plot(Mean_FD_Matrix(:,1),Mean_FD_Matrix(:,2))
                        title(['The Functioanl Mean of ERP with: ',NameBasisFunc,' Basis Function.'])
                        xlabel('Time (ms)') 
                        ylabel('Mean')
                        subplot(2,2,2)
                        plot(TimeIntervalSeleInput,mean(ADS_EData(SelectedChanIndex,:,:),3)')
                        title("The Mean of ERP")
                        xlabel('Time (ms)') 
                        ylabel('Mean')
                        subplot(2,2,3)
                        plot(STD_FD_Matrix(:,1),STD_FD_Matrix(:,2))
                        title(['The Functioanl Standard Deviation of ERP with: ',NameBasisFunc,' Basis Function.'])
                        xlabel('Time (ms)') 
                        ylabel('Standard Deviation')
                        subplot(2,2,4)
                        plot(TimeIntervalSeleInput,std(ADS_EData(SelectedChanIndex,:,:),0,3)')
                        title("The Standard Deviation of ERP")
                        xlabel('Time (ms)') 
                        ylabel('Standard Deviation')
                    end 
                       
                       if NSelectedChan > 1
                        fprintf('The number of channel is greater than 1.\n')  
                        fprintf('The functioal mean and standard deviation are plotted.\n')
                        figure()
                        subplot(2,2,1)
                        plot(Mean_FD_Matrix(:,1),Mean_FD_Matrix(:,2:(NSelectedChan+1)))
                        title(['The Functioanl Mean of ERP with: ',NameBasisFunc,' Basis Function.'])
                        xlabel('Time (ms)') 
                        ylabel('Mean')
                        subplot(2,2,2)
                        plot(TimeIntervalSeleInput,mean(ADS_EData(SelectedChanIndex,:,:),3)')
                        title("The Mean of ERP")
                        xlabel('Time (ms)') 
                        ylabel('Mean')
                        subplot(2,2,3)
                        plot(STD_FD_Matrix(:,1),STD_FD_Matrix(:,2:(NSelectedChan+1)))
                        title(['The Functioanl Standard Deviation of ERP with: ',NameBasisFunc,' Basis Function.'])
                        xlabel('Time (ms)') 
                        ylabel('Standard Deviation')
                        subplot(2,2,4)
                        plot(TimeIntervalSeleInput,std(ADS_EData(SelectedChanIndex,:,:),0,3)')
                        title("The Standard Deviation of ERP")
                        xlabel('Time (ms)') 
                        ylabel('Standard Deviation')
                       end 
                   end
      %%%%%%%%%%%%%%%% Calculate and Plot derivative from NODERIVATIVE  to 1 
                if Deri_PP_plot  == 1
                     fprintf('The derivative calculation is selectes.\n')  
                     SIZE_FD =  size(Mean_FD_Matrix);
                     Mean_FD_Matrix_Data = Mean_FD_Matrix(:,2:SIZE_FD(2));                     
                     Mean_FD_Matrix_Time = Mean_FD_Matrix(:,1);
      %%%%%%%%%%%%%%%% Basis Definition 
                if TypeBasisFunc == 1
      %%%%%%%%%%%%%%%%%%% B-Spline Basis   
                 wbasis = create_bspline_basis(rng, NumBasis, norder);
                 ADS_EData_Bspline_PP = double(squeeze(Mean_FD_Matrix_Data));
                 ADS_EData_Bspline_PP = ADS_EData_Bspline_PP';
                 ADS_EData_Bspline_PP = double(ADS_EData_Bspline_PP');
                 TimeIntervalSeleInput_BSpline_PP = double(Mean_FD_Matrix_Time');
                end 
                if TypeBasisFunc == 2
      %%%%%%%%%%%%%%%%%%% Fourier Basis   
                 wbasis = create_fourier_basis(rng, NumBasis);
                 ADS_EData_Bspline_PP = double(squeeze(Mean_FD_Matrix_Data));
                 ADS_EData_Bspline_PP = ADS_EData_Bspline_PP';
                 ADS_EData_Bspline_PP = double(ADS_EData_Bspline_PP');
                 TimeIntervalSeleInput_BSpline_PP = double(Mean_FD_Matrix_Time');
                end
      %%%%%%%%%%%%%%%% GCV for Lambda estimation 
                 if TickGCV == 1
                      fprintf('The smoothing of channel means for derivative with GCV start.\n')  
                           n = length(ADS_EData_Bspline_PP(1,:));
                      loglam = (-9:0.25:9)';
                        nlam = length(loglam);
                         Lfd = 4;
                    dfsave  = zeros(nlam,1);
                    gcvsave = zeros(nlam,n)';
                    MSEsave = zeros(nlam,n)';
                    for ilam=1:nlam
                        lambda = 10^loglam(ilam);
                        hgtfdPar = fdPar(wbasis, Lfd, lambda);
                        [hgtfd, df, gcv, coef,SSE] = smooth_basis(TimeIntervalSeleInput_BSpline_PP,ADS_EData_Bspline_PP,hgtfdPar);
                        accest = eval_fd(TimeIntervalSeleInput_BSpline_PP, hgtfd, 2);
                        dfsave(ilam)    = df;
                        gcvsave(:,ilam)   = gcv; 
                        MSEsave(:,ilam)   = SSE;
                    end
                    [MVALVE,MININDEX] =  min(mean(gcvsave));
                    LOGLAMMIN = loglam(MININDEX);
                    lambdaSelected  = 10^loglam(MININDEX);         % Minimum Lambda
                    hgtfdPar = fdPar(wbasis, Lfd, lambdaSelected);
                    ADS_EData_Bspline_PPP = smooth_basis(TimeIntervalSeleInput_BSpline_PP,ADS_EData_Bspline_PP,hgtfdPar);
                      fprintf('The smoothing of channel means for derivative with GCV end.\n')  

     %%%%%%%%%%%%%%%% Calculate and Plot derivative from NODERIVATIVE  to 1 
             if PPL_DER  == 1
                   fprintf('The derivative plot is selected.\n')  
                   Nder = NODERIVATIVE;
                   CM = jet(Nder);
                   for ii=1:Nder
                       fprintf(['The derivative plot ', num2str(ii) ,' is selected.\n'])  
                       Derive_FD =  deriv_fd(ADS_EData_Bspline_PPP,int2Lfd(ii));
                       curve_Deriv_FD = getcoef(Derive_FD);
                       curve_Deriv_FD_range = getbasisrange(getbasis(Derive_FD));
                       curve_Deriv_FD_nbasis = getnbasis(getbasis(Derive_FD));
                       curve_Deriv_FD_nx = max([501, 10*curve_Deriv_FD_nbasis+1]);
                       curve_Deriv_FD_x  = linspace(curve_Deriv_FD_range(1),curve_Deriv_FD_range(2),curve_Deriv_FD_nx)';
                       curve_Deriv_FD_Lfdobj = int2Lfd(int2Lfd(0));                           
                       curve_Deriv_FD_fdmat = eval_fd(curve_Deriv_FD_x,Derive_FD, curve_Deriv_FD_Lfdobj);
                       curve_Deriv_FD_max = max(curve_Deriv_FD_fdmat);
                       curve_Deriv_FD_min = min(curve_Deriv_FD_fdmat);
                       MinS = min(mean(curve_Deriv_FD_fdmat,2));
                       MaxS = max(mean(curve_Deriv_FD_fdmat,2));
                       MinT = min(min(TimeIntervalSeleInput),min(TimeIntervalSeleInput));
                       MaxT = max(max(TimeIntervalSeleInput),max(TimeIntervalSeleInput));    
                       figure()
                       plot(curve_Deriv_FD_x,mean(curve_Deriv_FD_fdmat,2),'color',CM(ii,:),'LineWidth',1.0,'LineStyle' ,'-')
                       title(['The Functional mean of ',num2str(ii),' derivative.'] )
                       line([MinT,MaxT], [0,0],'Color',[0 0 0]+0.05,'LineStyle','--')
                       axis([MinT-10 MaxT+10 (MinS) (MaxS)])
                       xlabel('Time (ms)') 
                       ylabel('Mean')
                   end                  
             end
                      
      %%%%%%%%%%%%%%%% Calculate and Plot Phase-Plane Plots                 
                   if PPL== 1
                       fprintf('The Phase-Plane Plot is selected.\n')  
                       Derive_FD_1 =  deriv_fd(ADS_EData_Bspline_PPP,int2Lfd(1));
                       curve_Deriv_FD_1 = getcoef(Derive_FD_1);
                       curve_Deriv_FD_range_1 = getbasisrange(getbasis(Derive_FD_1));
                       curve_Deriv_FD_nbasis_1 = getnbasis(getbasis(Derive_FD_1));
                       curve_Deriv_FD_nx_1 = max([501, 10*curve_Deriv_FD_nbasis_1+1]);
                       curve_Deriv_FD_x_1  = linspace(curve_Deriv_FD_range_1(1),curve_Deriv_FD_range_1(2),curve_Deriv_FD_nx_1)';
                       curve_Deriv_FD_Lfdobj_1 = int2Lfd(int2Lfd(0));                           
                       curve_Deriv_FD_fdmat_1 = eval_fd(curve_Deriv_FD_x_1,Derive_FD_1, curve_Deriv_FD_Lfdobj_1);
                       curve_Deriv_FD_max_1 = max(curve_Deriv_FD_fdmat_1);
                       curve_Deriv_FD_min_1 = min(curve_Deriv_FD_fdmat_1);
                       Derive_FD_2 =  deriv_fd(ADS_EData_Bspline_PPP,int2Lfd(2));
                       curve_Deriv_FD_2 = getcoef(Derive_FD_2);
                       curve_Deriv_FD_range_2 = getbasisrange(getbasis(Derive_FD_2));
                       curve_Deriv_FD_nbasis_2 = getnbasis(getbasis(Derive_FD_2));
                       curve_Deriv_FD_nx_2 = max([501, 10*curve_Deriv_FD_nbasis_2+1]);
                       curve_Deriv_FD_x_2  = linspace(curve_Deriv_FD_range_2(1),curve_Deriv_FD_range_2(2),curve_Deriv_FD_nx_2)';
                       curve_Deriv_FD_Lfdobj_2 = int2Lfd(int2Lfd(0));                           
                       curve_Deriv_FD_fdmat_2 = eval_fd(curve_Deriv_FD_x_2,Derive_FD_2, curve_Deriv_FD_Lfdobj_2);
                       curve_Deriv_FD_max_2 = max(curve_Deriv_FD_fdmat_2);
                       curve_Deriv_FD_min_2 = min(curve_Deriv_FD_fdmat_2);                                 
                       timeDomain = curve_Deriv_FD_x_1(:,1);
                       SizeTDomain = size(timeDomain);
                       indxTDomain = 1:30:SizeTDomain(1);
                       Xvar = mean(curve_Deriv_FD_fdmat_1,2);
                       Yvar = mean(curve_Deriv_FD_fdmat_2,2);
                       timeDomainC = num2str(round(timeDomain(indxTDomain),0));
                       ddd = size(timeDomainC);
                       indxTDomain_Zero = 1:1:sum(timeDomain < 0);
                       Xvar_Zero = mean(curve_Deriv_FD_fdmat_1,2);
                       Yvar_Zero = mean(curve_Deriv_FD_fdmat_2,2);
                       timeDomainC_Zero = num2str(round(timeDomain(indxTDomain_Zero),0));        
                       figure()
                       plot(Xvar,Yvar,'LineWidth',1.0,'LineStyle' ,'-','color','blue')
                       hold on
                       plot(Xvar_Zero(indxTDomain_Zero),Yvar_Zero(indxTDomain_Zero),'LineWidth',1.0,'LineStyle' ,'-','color','red')
                       hold on
                       plot(Xvar(indxTDomain),Yvar(indxTDomain),'*','color','black')
                       hold on 
                       text(Xvar(indxTDomain),Yvar(indxTDomain),timeDomainC(:,1:ddd(2)),'FontSize', 10)
                       hold on 
                       plot(Xvar(1),Yvar(1),'*','color','black','LineWidth',10.0)
                       hold on 
                       plot(Xvar(SizeTDomain(1)),Yvar(SizeTDomain(1)),'*','color','black','LineWidth',10.0)
                       xlabel('\fontsize{13} Velocity (\mu. volt/ms)')
                       ylabel('\fontsize{13} Acceleration (\mu. volt/msˆ2)')
                       title(['The Phase-Plane Plot'] )
                  end                  
            end
        end 
     
      if FPCA_STATUS == 1
      %%%%%%%%%%%%%%%% Preprocessing Steps
         fprintf('The Functioanl Principal Component Analysis is selected.\n')  
         SIZE_FD_FPCA =  size(Mean_FD_Matrix);
         Mean_FD_Matrix_Data_FDA = Mean_FD_Matrix(:,2:SIZE_FD_FPCA(2));                     
         Mean_FD_Matrix_Time_FDA = Mean_FD_Matrix(:,1);
      %%%%%%%%%%%%%%%% Basis Definition 
                if TypeBasisFunc == 1
      %%%%%%%%%%%%%%%%%%% B-Spline Basis   
                 wbasis_FCPA = create_bspline_basis(rng, NumBasis, norder);
                 ADS_EData_Bspline_FPCA = double(squeeze(Mean_FD_Matrix_Data_FDA));
                 ADS_EData_Bspline_FPCA = ADS_EData_Bspline_FPCA';
                 ADS_EData_Bspline_FPCA = double(ADS_EData_Bspline_FPCA');
                 TimeIntervalSeleInput_BSpline_FPCA = double(Mean_FD_Matrix_Time_FDA');
                end 
                if TypeBasisFunc == 2
      %%%%%%%%%%%%%%%%%%% Fourier Basis   
                 wbasis_FCPA = create_fourier_basis(rng, NumBasis);
                 ADS_EData_Bspline_FPCA = double(squeeze(Mean_FD_Matrix_Data_FDA));
                 ADS_EData_Bspline_FPCA = ADS_EData_Bspline_FPCA';
                 ADS_EData_Bspline_FPCA = double(ADS_EData_Bspline_FPCA');
                 TimeIntervalSeleInput_BSpline_FPCA = double(Mean_FD_Matrix_Time_FDA');
                end
      %%%%%%%%%%%%%%%% GCV for Lambda estimation 
                 if TickGCV == 1
                      fprintf('The smoothing of channel means for FPCA with GCV start.\n')  
                           n = length(ADS_EData_Bspline_FPCA(1,:));
                      loglam = (-9:0.25:9)';
                        nlam = length(loglam);
                         Lfd = 4;
                    dfsave  = zeros(nlam,1);
                    gcvsave = zeros(nlam,n)';
                    MSEsave = zeros(nlam,n)';
                    for ilam=1:nlam
                        lambda = 10^loglam(ilam);
                        hgtfdPar = fdPar(wbasis_FCPA, Lfd, lambda);
                        [hgtfd, df, gcv, coef,SSE] = smooth_basis(TimeIntervalSeleInput_BSpline_FPCA,ADS_EData_Bspline_FPCA,hgtfdPar);
                        accest = eval_fd(TimeIntervalSeleInput_BSpline_FPCA, hgtfd, 2);
                        dfsave(ilam)    = df;
                        gcvsave(:,ilam)   = gcv; 
                        MSEsave(:,ilam)   = SSE;
                    end
                    [MVALVE,MININDEX] =  min(mean(gcvsave));
                    LOGLAMMIN = loglam(MININDEX);
                    lambdaSelected  = 10^loglam(MININDEX);         % Minimum Lambda
                    hgtfdPar = fdPar(wbasis, Lfd, lambdaSelected);
                    ADS_EData_Bspline_FPCA = smooth_basis(TimeIntervalSeleInput_BSpline_FPCA,ADS_EData_Bspline_FPCA,hgtfdPar);
                    fprintf('The smoothing of channel means for FPCA with GCV end.\n')  
                 
      %%%%%%%%%%%%%%%% Plot Functional PCA                  
                   nharm  = FPCA_nFPCA;
                   EEG_FPCA = pca_fd(ADS_EData_Bspline_FPCA, nharm);
                   fprintf(['The number of selected FPCA are ',  num2str(nharm) ,' .\n'])  
                   rpt =   [repmat(SelectedChnInd,nharm,1),(1:nharm)',EEG_FPCA.varprop*100,cumsum(EEG_FPCA.varprop*100)];   
                   tableReport = array2table(rpt,'VariableNames',{'Channel','Number_FPCAs','FVE','Cumulative_FVE'})
          
                   
      %%%%%%%%%%%%%%%% FPCA Eigenfunction Plot                
                         if FPCA_PlotFPCA == 1
                               CM = jet(nharm);
                                for iii=1:nharm
                                       fprintf(['The functioanl eigenfunction ',  num2str(iii) ,' is plotted .\n'])  
                                       varExp = EEG_FPCA.varprop;
                                       harmfd_iii = EEG_FPCA.harmfd(iii);
                                       curve_iii = getcoef(harmfd_iii);
                                       curve_iii_range = getbasisrange(getbasis(harmfd_iii));
                                       curve_iii_nbasis   = getnbasis(getbasis(harmfd_iii));
                                       curve_iii_nx = max([501, 10*curve_iii_nbasis+1]);
                                       curve_iii_x     = linspace(curve_iii_range(1),curve_iii_range(2),curve_iii_nx)';
                                       curve_iii_Lfdobj = int2Lfd(int2Lfd(0));                           
                                       curve_iii_fdmat = eval_fd(curve_iii_x,harmfd_iii, curve_iii_Lfdobj);
                                       max_iii = max(curve_iii);
                                       min_iii = min(curve_iii);
                                       MinS = min(min_iii);
                                       MaxS = max(max_iii);
                                       MinT = min(min(TimeIntervalSeleInput),min(TimeIntervalSeleInput));
                                       MaxT = max(max(TimeIntervalSeleInput),max(TimeIntervalSeleInput));
                                       figure()
                                       plot(curve_iii_x,curve_iii_fdmat,'color',CM(iii,:),'LineWidth',1.5,'LineStyle' ,'-')
                                       title(['The PCA function ', num2str(iii) ,' (Percentage of variability ', num2str(round(varExp(iii)*100,2)) ,' %)'] )
                                       line([MinT,MaxT], [0,0],'Color',[0 0 0]+0.05,'LineStyle','--')
                                       axis([MinT-10 MaxT+10 MinS-0.025 MaxS+0.025])
                                       legend({[num2str(iii) ,' FPCA ']})
                                       xlabel('Time (ms)') 
                                       ylabel('Harmonic')
                                   end 
                         end
      %%%%%%%%%%%%%%%% FPCA Scree Plot               
                       if FPCA_PlotFPCA_SCREE == 1
                          fprintf(['The scree plot is selected. \n'])  
                              figure()
                              plot(EEG_FPCA.varprop,'-O')
                              xlabel('Eigenvalue Number')
                              ylabel('Fraction of Variance Explained %')
                              title(['Total FVE  (' num2str(round(sum(EEG_FPCA.varprop)*100,2)) '% with ', num2str(nharm) , ' FPCAs)'])
                       end
                       
     %%%%%%%%%%%%%%%% Name of Channels              
                       n_Chans = size(SelectedChanIndex);
                       strName_label = cell(n_Chans(2),1);
                       for hh = 1:n_Chans(2)
                          strName_label{hh,1} = string(SelectedChan(hh).labels);
                       end 
                       
     %%%%%%%%%%%%%%%% FPCA Scores Plot             
                       if FPCA_Scores == 1
                         fprintf(['The fpca scores plot is selected. \n'])  
                          varExp = EEG_FPCA.varprop;
                                  for ii = 1:FPCA_nFPCA
                                      for jj = 1:FPCA_nFPCA
                                          if ii < jj                            
                                             figure()
                                             fprintf(['The fpca scores plot of ' , num2str(ii) , ' against ', num2str(jj)  ,' . \n'])  
                                              plot(EEG_FPCA.harmscr(:,ii),EEG_FPCA.harmscr(:,jj),".","MarkerSize",14,'color','black')                                           
                                              line([0, 0], [min(EEG_FPCA.harmscr(:,jj))-100,max(EEG_FPCA.harmscr(:,jj))+100],'Color','red','LineStyle','--')
                                              line([min(EEG_FPCA.harmscr(:,ii))-100,max(EEG_FPCA.harmscr(:,ii))+100], [0,0],'Color','red','LineStyle','--')
                                              title(['The ', num2str(ii) ,' against ' ,num2str(jj), ' Harmonic scores function of ', SelectedChanLabels ])
                                              xlabel(['Scores ', num2str(ii), ' of ',SelectedChanLabels, '( ',num2str(round(varExp(ii),3)*100),' %)']) 
                                              ylabel(['Scores ', num2str(jj), ' of ',SelectedChanLabels, '( ',num2str(round(varExp(jj),3)*100),' %)'])
                                              axis([ (min(EEG_FPCA.harmscr(:,ii))-20) (max(EEG_FPCA.harmscr(:,ii))+20) (min(EEG_FPCA.harmscr(:,jj))-20) (max(EEG_FPCA.harmscr(:,jj))+20) ])
                                             text(EEG_FPCA.harmscr(:,ii),EEG_FPCA.harmscr(:,jj),strName_label,'Color', 'blue','FontSize',10) 
                                          end 
                                        end
                                  end
                       end  
      %%%%%%%%%%%%%%%% FPCA Scores Plot             
                         if FPCA_Scores_3D == 1     
                                              figure();
                                              fprintf('The fpca 3D plot is selected. \n')  
                                              scatter3(EEG_FPCA.harmscr(:,1),EEG_FPCA.harmscr(:,2),EEG_FPCA.harmscr(:,3),'MarkerFaceColor',[0 .75 .75])
                                              title('The frist three FPCA Scores')
                                              text(EEG_FPCA.harmscr(:,1),EEG_FPCA.harmscr(:,2),EEG_FPCA.harmscr(:,3),strName_label,'Color', 'red','FontSize',10) 
                                              xlabel(['Scores FPCA ', num2str(1), '( ',num2str(round(varExp(1),3)*100),' %)']) 
                                              ylabel(['Scores FPCA ', num2str(2), '( ',num2str(round(varExp(2),3)*100),' %)'])
                                              zlabel(['Scores FPCA ', num2str(3), '( ',num2str(round(varExp(3),3)*100),' %)'])
                         end
       %%%%%%%%%%%%%%%% FPCA VaRIMAX Plot        
                          if FPCA_PlotFPCA_VARMAX == 1
                            fprintf('The VARIMAX rotation is selected.')  
                              EEG_FPCA_VARMAX = varmx_pca(EEG_FPCA);
                               CM = jet(nharm);
                                for iii=1:nharm
                                       fprintf(['The VARIMAX rotation of functioanl eigenfunction ',  num2str(iii) ,' is plotted .\n'])  
                                       varExp = EEG_FPCA_VARMAX.varprop;
                                       harmfd_iii = EEG_FPCA_VARMAX.harmfd(iii);
                                       curve_iii = getcoef(harmfd_iii);
                                       curve_iii_range = getbasisrange(getbasis(harmfd_iii));
                                       curve_iii_nbasis   = getnbasis(getbasis(harmfd_iii));
                                       curve_iii_nx = max([501, 10*curve_iii_nbasis+1]);
                                       curve_iii_x     = linspace(curve_iii_range(1),curve_iii_range(2),curve_iii_nx)';
                                       curve_iii_Lfdobj = int2Lfd(int2Lfd(0));                           
                                       curve_iii_fdmat = eval_fd(curve_iii_x,harmfd_iii, curve_iii_Lfdobj);
                                       max_iii = max(curve_iii);
                                       min_iii = min(curve_iii);
                                       MinS = min(min_iii);
                                       MaxS = max(max_iii);
                                       MinT = min(min(TimeIntervalSeleInput),min(TimeIntervalSeleInput));
                                       MaxT = max(max(TimeIntervalSeleInput),max(TimeIntervalSeleInput));
                                       figure()
                                       plot(curve_iii_x,curve_iii_fdmat,'color',CM(iii,:),'LineWidth',1.5,'LineStyle' ,'-')
                                       title(['The Varimax Rotation of PCA function ', num2str(iii)  ,' (Percentage of variability ', num2str(round(varExp(iii)*100,2)) ,' %)'] )
                                       line([MinT,MaxT], [0,0],'Color',[0 0 0]+0.05,'LineStyle','--')
                                       axis([MinT-10 MaxT+10 MinS-0.025 MaxS+0.025])
                                       legend({[num2str(iii) ,' FPCA ']})
                                       xlabel('Time (ms)') 
                                       ylabel('Harmonic')
                                end
                          end
         %%%%%%%%%%%%%%%% FPCA Varimax Scree Plot     
                              if FPCA_VARIMAX_PlotFPCA_SCREE == 1
                                    fprintf(['The scree plot for VARIMAX Rotation is selected. \n'])  
                                    figure()
                                    plot(EEG_FPCA_VARMAX.varprop,'-O')
                                    xlabel('Eigenvalue Number')
                                    ylabel('Fraction of Variance Explained %')
                                    title(['Total FVE in Varimax (' num2str(round(sum(EEG_FPCA_VARMAX.varprop)*100,2)) ' with ', num2str(nharm) , ' FPCAs)'])
                              end
          %%%%%%%%%%%%%%%% FPCA Varimax Scores Plot     
                             if FPCA_Scores_VARMAX == 1 
                                     fprintf(['The VARIMAX rotation of fpca scores plot is selected. \n'])  
                                     varExp_VARMAX = EEG_FPCA_VARMAX.varprop;
                                              for ii = 1:FPCA_nFPCA
                                                  for jj = 1:FPCA_nFPCA
                                                     if ii < jj  
                                                        fprintf(['The VAIMAX rotation of fpca scores plot of ' , num2str(ii) , ' against ', num2str(jj)  ,' . \n'])  
                                                        figure()
                                                        plot(EEG_FPCA_VARMAX.harmscr(:,ii),EEG_FPCA_VARMAX.harmscr(:,jj),".","MarkerSize",14,'color','black')                                           
                                                        line([0, 0], [min(EEG_FPCA_VARMAX.harmscr(:,jj))-100,max(EEG_FPCA_VARMAX.harmscr(:,jj))+100],'Color','red','LineStyle','--')
                                                        line([min(EEG_FPCA_VARMAX.harmscr(:,ii))-100,max(EEG_FPCA_VARMAX.harmscr(:,ii))+100], [0,0],'Color','red','LineStyle','--')
                                                        title(['The ', num2str(ii) ,' against ' ,num2str(jj), ' Varimax Harmonic scores function' ])
                                                        xlabel(['Scores ', num2str(ii), ' of ',SelectedChanLabels, '( ',num2str(round(varExp(ii),3)*100),' %)']) 
                                                        ylabel(['Scores ', num2str(jj), ' of ',SelectedChanLabels, '( ',num2str(round(varExp(jj),3)*100),' %)'])
                                                        axis([ (min(EEG_FPCA_VARMAX.harmscr(:,ii))-20) (max(EEG_FPCA_VARMAX.harmscr(:,ii))+20) (min(EEG_FPCA_VARMAX.harmscr(:,jj))-20) (max(EEG_FPCA_VARMAX.harmscr(:,jj))+20) ])
                                                        text(EEG_FPCA_VARMAX.harmscr(:,ii),EEG_FPCA_VARMAX.harmscr(:,jj),strName_label,'Color', 'blue','FontSize',8) 
                                                        end 
                                                     end
                                              end
                             end
                                              
            %%%%%%%%%%%%%%%% FPCA Scores Plot             
                         if FPCA_Scores_VARMAX_3D == 1    
                              fprintf('The fpca 3D plot is selected. \n')  
                               figure();
                               scatter3(EEG_FPCA_VARMAX.harmscr(:,1),EEG_FPCA_VARMAX.harmscr(:,2),EEG_FPCA_VARMAX.harmscr(:,3),'MarkerFaceColor',[0 .75 .75])
                               title('The Frist three Varimax FPCA Scores')
                               text(EEG_FPCA_VARMAX.harmscr(:,1),EEG_FPCA_VARMAX.harmscr(:,2),EEG_FPCA_VARMAX.harmscr(:,3),strName_label,'Color', 'red','FontSize',10) 
                               xlabel(['Scores FPCA ', num2str(1), '( ',num2str(round(varExp_VARMAX(1),3)*100),' %)']) 
                               ylabel(['Scores FPCA ', num2str(2), '( ',num2str(round(varExp_VARMAX(2),3)*100),' %)'])
                               zlabel(['Scores FPCA ', num2str(3), '( ',num2str(round(varExp_VARMAX(3),3)*100),' %)'])
                              end   
         end             
end