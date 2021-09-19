% pop_fda_2 - perform epoch analysis with functional principcal component  
%             analsysi (FPCA). It also smooth curves with B-Spline and
%             Fourier basis functions. The estimation of penalty term is
%             obtained with generalized cross validavtion (GCV). 
% Usage:
%  pop_fda_2(EEG); % pop up window asking users to select method
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
%   pop_fda_2(EEG); % pop up window asking users to select method
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


function com = pop_fda_2(EEG, varargin)
%%%% Tools -> EPOCH Analysis -> FPCA

fprintf('EPOCH Analysis \n');
fprintf('Functional Principal Component Analysis (FPCA) \n');

com = '';
if nargin < 1
    help pop_sourcereconstruction;
    return
end

if nargin < 2
    %%%% GUI 
    g = [1 0.5 0.5];
    geometry = { g g g g g g g g g g g g g g g g};
    uilist = { ...
      { 'Style', 'text', 'string', 'Channel Selection', 'fontweight', 'bold' } {} {} ...
      { 'Style', 'text', 'string', '   Channel number', 'FontAngle', 'italic'  } ...
      { 'Style', 'edit', 'string', '1:32' 'tag' 'SelectedChanIndex'} {} ...
      ...
      { 'Style', 'text', 'string', '   Time limits [min max] (msec)', 'FontAngle', 'italic' } ...
      { 'Style', 'edit', 'string', '-0.2 0.800 ' 'tag' 'TIS' } {} ...
      ...
      { 'Style', 'text', 'string', 'Preprocessing (Step 1) ', 'fontweight', 'bold' } {} {} ...
      { 'Style', 'text', 'string', '   Choose Basis Function (Default: B-Spline)', 'FontAngle', 'italic' } ...
      { 'Style', 'popupmenu', 'string', 'B-Spline|Fourier' , 'tag' 'TypeBasisFunc' } {} ...
      { 'Style', 'text', 'string', '   Number of Basis', 'FontAngle', 'italic' } ...
      { 'Style', 'edit', 'string', '120' 'tag' 'NB'  } {} ...
      ...
      { 'Style', 'text', 'string', '   Order of Basis','FontAngle', 'italic' } ...
      { 'Style', 'edit', 'string', '6' 'tag' 'NORD' } {} ...
      { 'Style', 'text', 'string', '   Parameter Estimation', 'FontAngle', 'italic' } ...
      { 'Style', 'checkbox', 'string' 'GCV' 'value' 1 'tag' 'TickGCV' } ...
      { 'Style', 'checkbox', 'string' 'Plot GCV' 'value' 1 'tag' 'PlotGCV' } ...
      { 'Style', 'text', 'string', 'Functional Principal Components Analysis (FPCA)', 'fontweight', 'bold' } {} {} ...
      { 'Style', 'text', 'string', '   FPCA Status',  'FontAngle', 'italic' } ...
      { 'Style', 'checkbox', 'string' 'Done' 'value' 1 'tag' 'FPCA_STATUS' } {} ...
      { 'Style', 'text', 'string', '   Number of FPCA', 'FontAngle', 'italic' } ...
      { 'Style', 'edit', 'string', '4' 'tag' 'FPCA_nFPCA'} {} ...
      { 'Style', 'text', 'string', '   FPCA Descriptive', 'FontAngle', 'italic' } ...
      { 'Style', 'checkbox', 'string' 'Plot Scree plot' 'value' 1 'tag' 'FPCA_PlotFPCA_SCREE' } ...
      { 'Style', 'checkbox', 'string' 'Plot Eigefunction' 'value' 0 'tag' 'FPCA_PlotFPCA' } ...
      { 'Style', 'text', 'string', '   Rotations', 'FontAngle', 'italic' } ...
      { 'Style', 'checkbox', 'string' 'Varimax Rotation' 'value' 1 'tag' 'FPCA_PlotFPCA_VARMAX' }  ...
      { 'Style', 'checkbox', 'string' 'Plot Scree plot' 'value' 0 'tag' 'FPCA_VARIMAX_PlotFPCA_SCREE' } ...      
      { 'Style', 'text', 'string', '   Plot Scores', 'FontAngle', 'italic' } ...
      { 'Style', 'checkbox', 'string' 'Scores' 'value' 0 'tag' 'FPCA_Scores' } ...
      { 'Style', 'checkbox', 'string' 'Varimax Scores' 'value' 0 'tag' 'FPCA_Scores_VARMAX'} ...   
      { 'Style', 'text', 'string', 'General Setting', 'fontweight', 'bold' } {} {} ...
      { 'Style', 'text', 'string', '   Plot Setting', 'FontAngle', 'italic' } ...
      { 'Style', 'checkbox', 'string' 'Plot All Smooth Curve' 'value' 0 'tag' 'AllSmoothPlot' } {} ...
      };

      %%%% Calling (GUI)
     [ tmp1 tmp2 strhalt structout ] = inputgui( geometry, uilist, ...
           'pophelp(''pop_newtimef'');', 'EPOCH Analysis | Functional Principal Component Analysis for EEGLAB - Version 0.1');

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
            SelectedChanIndex = str2num(structout.SelectedChanIndex);
            SelectedChan = ChanLocs(:,SelectedChanIndex) ;
            NSelectedChan = length(SelectedChanIndex);

            %%% Event Selection
            % SelectedEventIndex = M;
            % SelectedEvent = Events(:,SelectedEventIndex);

            %%% Time Interval Selection 
            TIS = str2num(structout.TIS) ;
            TimeIntervalSele = TIS; %%% Time Interval Selection (defualt: -200, 800 ms)

            TimeLO = min(TimeIntervalSele) * 1.0e+03 ;
            TimeUP = max(TimeIntervalSele) * 1.0e+03 ;

            TimeIntervalSeleInput      = Times(Times(1,:) >= TimeLO & Times(1,:) <= TimeUP);
            TimeIntervalSeleIndexInput = Times(1,:) >= TimeLO & Times(1,:) <= TimeUP;

            %%% Analysis Data Set with selected time intervals
            ADS_EData = EData(:,TimeIntervalSeleIndexInput(1,:),:);

            %%% Analysis 

            %%%% Preprocessing
            %%% B-Spline
            %% Input
            NumBasis = str2num(structout.NB); % Number of Splines
            norder   = str2num(structout.NORD); % Number of Order
            TickGCV  = structout.TickGCV; % If GCV is ticked (TRUE-Default), the parameters are optimised with GCV, O.W. it is FALSE.
            PlotGCV = structout.PlotGCV; % Plot GCV against Lambda Status (Deafult=0), If PlotGCV=1, it plots
            rng = [min(TimeIntervalSeleInput),max(TimeIntervalSeleInput)];  % Range of B-Spline
            ADS_EData_Bspline_1_out=cell(NSelectedChan,1,nTrials);

            FPCA_STATUS = structout.FPCA_STATUS; % Doing FPCA (Default = 0), if FPCA_STATUS=1, plot FPCA. 
            FPCA_nFPCA = str2num(structout.FPCA_nFPCA) ; % Number of FPCA
            FPCA_PlotFPCA = structout.FPCA_PlotFPCA; % Plot Eigenfunctions (Deafult = 0), if FPCA_PlotFPCA=1, Plot EigenFunctions
            FPCA_PlotFPCA_VARMAX =structout.FPCA_PlotFPCA_VARMAX; % Plot Varimax Eigenfunctions (Deafult = 0), if FPCA_PlotFPCA=1, Plot Varimax EigenFunctions
            FPCA_PlotFPCA_SCREE = structout.FPCA_PlotFPCA_SCREE; % Plot Scree plot (Deafult = 0), if FPCA_PlotFPCA_SCREE=1, Plot Scree Plot
            FPCA_VARIMAX_PlotFPCA_SCREE = structout.FPCA_VARIMAX_PlotFPCA_SCREE;
            FPCA_Scores = structout.FPCA_Scores;
            FPCA_Scores_VARMAX = structout.FPCA_Scores_VARMAX;
            TypeBasisFunc = structout.TypeBasisFunc;
            
            
            AllSmoothPlot = structout.AllSmoothPlot; % Plot all smoothed curves (default = 0 ) , if AllSmoothPlot=1, plot all smoothed curves. 
            %% Anlaysis Code
            for nch = 1:NSelectedChan
                SelectedChnInd = SelectedChanIndex(nch);
                SelectedChanLabels = SelectedChan(nch).labels;
                ADS_EData_Bspline = ADS_EData(SelectedChnInd,:,:);

                 if TypeBasisFunc == 1
                 %%% Defining B-Spline Basis  
                 %%% Preparing Data
                 wbasis = create_bspline_basis(rng, NumBasis, norder);
                 ADS_EData_Bspline = double(squeeze(ADS_EData_Bspline));
                 ADS_EData_Bspline = ADS_EData_Bspline';
                 ADS_EData_Bspline = double(ADS_EData_Bspline');
                 TimeIntervalSeleInput_BSpline = double(TimeIntervalSeleInput');
                end 
                if TypeBasisFunc == 2
                 %%% Defining Fourier Basis  
                 %%% Preparing Data
                 wbasis = create_fourier_basis(rng, NumBasis);
                 ADS_EData_Bspline = double(squeeze(ADS_EData_Bspline));
                 ADS_EData_Bspline = ADS_EData_Bspline';
                 ADS_EData_Bspline = double(ADS_EData_Bspline');
                 TimeIntervalSeleInput_BSpline = double(TimeIntervalSeleInput');
                end
            

                 if TickGCV == 1
                           n = length(ADS_EData_Bspline(1,:));
                      loglam = (-10.5:0.25:10.5)'; %  set up the range of log lambda values
                        nlam = length(loglam);
                         Lfd = 4;

                    dfsave  = zeros(nlam,1);
                    gcvsave = zeros(nlam,n)';
                    MSEsave = zeros(nlam,n)';

                    %  loop through the log lambda values (Modify Later this part for parameters)
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

                    lambdaSelected  = 10^loglam(MININDEX);
                    hgtfdPar = fdPar(wbasis, Lfd, lambdaSelected);
                    ADS_EData_Bspline_1 = smooth_basis(TimeIntervalSeleInput_BSpline,ADS_EData_Bspline,hgtfdPar);
                    fprintf(['*** Channel Number: ', num2str(SelectedChnInd),' , Minimum log Lambda: ',num2str(LOGLAMMIN), ' , Minimum GCV: ', num2str(MVALVE)]);

                    if PlotGCV == 1
                    %%% Plot GCV (Modify)
                   figure()
                    plot(loglam,mean(gcvsave)' )        
                    axis([ (min(loglam)-3) (max(loglam)+3) (min(mean(gcvsave)-5)) (max(mean(gcvsave)+5)) ])
                    suptitle(['Channel Name: ',SelectedChanLabels])
                    title(['The minimum GCV (', num2str(MVALVE) ,') at Log Lambda:' num2str(LOGLAMMIN) ])
                    yy = MVALVE;
                    xx = loglam(MININDEX);
                   line([xx, xx], [min(mean(gcvsave))-20 ,max(mean(gcvsave))+20],'Color','red','LineStyle','--')
                   line([min(loglam)-10, max(loglam)+10], [yy, yy],'Color','red','LineStyle','--')
                   ylabel('GCV') 
                   xlabel('Logarithm Lambda') 
                   end

                    %%% Plot MSE
                    %%%% plot(loglam,mean(MSEsave)) --> Based on Values
                    %%%% plot(mean(MSEsave))  %% --> Based on Index

                 elseif TickGCV == 0
                    %%% Smoothing Data
                    [ADS_EData_Bspline_1,df_1, gcv_1,coef_1,SSE_1] = smooth_basis(TimeIntervalSeleInput_BSpline,ADS_EData_Bspline,wbasis);
                     LOGLAMMIN =0;
                 end

                    %%% Output
                    EEGDatafd_names{1} = 'Time (ms)';
                    EEGDatafd_names{2} = SelectedChanLabels ;
                    EEGDatafd_names{3} = '\mu. volt';
                    ADS_EData_Bspline_1 = putnames(ADS_EData_Bspline_1, EEGDatafd_names);

                   %%%% Plot each Curve Separatley  
                   %%% plotfit_fd(ADS_EData_Bspline, TimeIntervalSeleInput_BSpline, ADS_EData_Bspline_1)
                   if AllSmoothPlot == 1
                   AllSmoothPlot_1 = getcoef(ADS_EData_Bspline_1);
                   AllSmoothPlot_range = getbasisrange(getbasis(ADS_EData_Bspline_1));
                   AllSmoothPlot_nbasis = getnbasis(getbasis(ADS_EData_Bspline_1));
                   AllSmoothPlot_nx = max([501, 10*AllSmoothPlot_nbasis+1]);
                   AllSmoothPlot_x  = linspace(AllSmoothPlot_range(1),AllSmoothPlot_range(2),AllSmoothPlot_nx)';
                   AllSmoothPlot_Lfdobj = int2Lfd(int2Lfd(0));                           
                   AllSmoothPlot_fdmat = eval_fd(AllSmoothPlot_x,ADS_EData_Bspline_1, AllSmoothPlot_Lfdobj);
                   AllSmoothPlot_max = max(AllSmoothPlot_fdmat);
                   AllSmoothPlot_min = min(AllSmoothPlot_fdmat);
                   
                   MaxS = max(AllSmoothPlot_max);
                   MinS = min(AllSmoothPlot_min);
                   
                   if TypeBasisFunc == 1 
                    MName = 'B-Spline';
                   end
                   
                   if TypeBasisFunc == 2 
                    MName = 'Fourier';
                   end

                   MinT = min(TimeIntervalSeleInput_BSpline);
                   MaxT = max(TimeIntervalSeleInput_BSpline);
                   
                   figure();
                   plot(AllSmoothPlot_x ,AllSmoothPlot_fdmat)
                   title(['The  Smoothing Functions of ', SelectedChanLabels ,' is:  ', MName ] )
                   line([MinT,MaxT], [0,0],'Color',[0 0 0]+0.05,'LineStyle','--')
                   suptitle(['Channel Name: ',SelectedChanLabels,' (Best Log Lambda:', num2str(LOGLAMMIN),')'])
                   axis([MinT-10 MaxT+10 MinS-25 MaxS+25])
                   ylabel('{\mu} Volt')
                   xlabel('Time (ms)')
                   end
                   %%% Output
                   %% Save all fd Object 
                   ADS_EData_Bspline_1_out(nch,1,:) = fd2cell(ADS_EData_Bspline_1);
                

                   if FPCA_STATUS == 1
                   %% Plot Functional PCA   
                   nharm  = FPCA_nFPCA;
                   EEG_FPCA = pca_fd(ADS_EData_Bspline_1, nharm);
                   rpt =   [repmat(SelectedChnInd,nharm,1),(1:nharm)',EEG_FPCA.varprop*100,cumsum(EEG_FPCA.varprop*100)];
                  
                   tableReport = array2table(rpt,'VariableNames',{'Channel','Number_FPCAs','FVE','Cumulative_FVE'})
          
                         %%% FPCA Eigenfunction Plot
                         if FPCA_PlotFPCA == 1
                               CM = jet(nharm);
                                for iii=1:nharm
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
                                       title(['The PCA function ', num2str(iii) ,' of ', SelectedChanLabels ,' (Percentage of variability ', num2str(round(varExp(iii)*100,2)) ,' %)'] )
                                       line([MinT,MaxT], [0,0],'Color',[0 0 0]+0.05,'LineStyle','--')
                                       axis([MinT-10 MaxT+10 MinS-0.025 MaxS+0.025])
                                       legend({[num2str(iii) ,' FPCA ', SelectedChanLabels]})
                                       xlabel('Time (ms)') 
                                       ylabel('Harmonic')
                                       % plot(EEG_FPCA.harmfd(iii))
                                       % plot_pca_fd(EEG_FPCA)
                                   end 
                         end

                    %%% FPCA Scree Plot
                       if FPCA_PlotFPCA_SCREE == 1
                            %  plot Scree plot
                              figure()
                              plot(EEG_FPCA.varprop,'-O')
                              xlabel('Eigenvalue Number')
                              ylabel('Fraction of Variance Explained %')
                              suptitle(['Channel Name: ',SelectedChanLabels])
                              title(['Total FVE  (' num2str(round(sum(EEG_FPCA.varprop)*100,2)) ' with ', num2str(nharm) , ' FPCAs)'])
                       end
                       
                  %%% FPCA Scores Plot
                       if FPCA_Scores == 1
                          varExp = EEG_FPCA.varprop;
                                  for ii = 1:FPCA_nFPCA
                                      for jj = 1:FPCA_nFPCA
                                          if ii < jj                            
                                             figure()
                                             strName = num2str((1:size(EEG_FPCA.harmscr(:,ii)))');
                                              plot(EEG_FPCA.harmscr(:,ii),EEG_FPCA.harmscr(:,jj),".","MarkerSize",14,'color','black')                                           
                                              line([0, 0], [min(EEG_FPCA.harmscr(:,jj))-100,max(EEG_FPCA.harmscr(:,jj))+100],'Color','red','LineStyle','--')
                                              line([min(EEG_FPCA.harmscr(:,ii))-100,max(EEG_FPCA.harmscr(:,ii))+100], [0,0],'Color','red','LineStyle','--')
                                              title(['The ', num2str(ii) ,' against ' ,num2str(jj), ' Harmonic scores function of ', SelectedChanLabels ])
                                              xlabel(['Scores ', num2str(ii), ' of ',SelectedChanLabels, '( ',num2str(round(varExp(ii),3)*100),' %)']) 
                                              ylabel(['Scores ', num2str(jj), ' of ',SelectedChanLabels, '( ',num2str(round(varExp(jj),3)*100),' %)'])
                                              axis([ (min(EEG_FPCA.harmscr(:,ii))-100) (max(EEG_FPCA.harmscr(:,ii))+100) (min(EEG_FPCA.harmscr(:,jj))-100) (max(EEG_FPCA.harmscr(:,jj))+100) ])
                                             text(EEG_FPCA.harmscr(:,ii),EEG_FPCA.harmscr(:,jj),strName,'Color', 'blue','FontSize',8) 
                                          end 
                                        end
                                  end
                       end  
                           
                           
                  %%% FPCA VaRIMAX Plot
                          if FPCA_PlotFPCA_VARMAX == 1
                              %  Varimax rotation
                              EEG_FPCA_VARMAX = varmx_pca(EEG_FPCA);
                               CM = jet(nharm);
                                for iii=1:nharm
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
                                       title(['The Varimax Rotation of PCA function ', num2str(iii) ,' of ', SelectedChanLabels ,' (Percentage of variability ', num2str(round(varExp(iii)*100,2)) ,' %)'] )
                                       line([MinT,MaxT], [0,0],'Color',[0 0 0]+0.05,'LineStyle','--')
                                       axis([MinT-10 MaxT+10 MinS-0.025 MaxS+0.025])
                                       legend({[num2str(iii) ,' FPCA ', SelectedChanLabels]})
                                       xlabel('Time (ms)') 
                                       ylabel('Harmonic')
                                end
                              
                          %%% FPCA Varimax Scree Plot
                              if FPCA_VARIMAX_PlotFPCA_SCREE == 1
                                 %  plot Scree plot
                                    figure()
                                    plot(EEG_FPCA_VARMAX.varprop,'-O')
                                    xlabel('Eigenvalue Number')
                                    ylabel('Fraction of Variance Explained %')
                                    suptitle(['Channel Name: ',SelectedChanLabels])
                                    title(['Total FVE in Varimax (' num2str(round(sum(EEG_FPCA_VARMAX.varprop)*100,2)) ' with ', num2str(nharm) , ' FPCAs)'])
                              end
                          
                  %%% FPCA Varimax Scores Plot
                              if FPCA_Scores_VARMAX == 1 
                                     varExp = EEG_FPCA_VARMAX.varprop;
                                              for ii = 1:FPCA_nFPCA
                                                  for jj = 1:FPCA_nFPCA
                                                     if ii < jj                            
                                                        figure()
                                                        strName = num2str((1:size(EEG_FPCA_VARMAX.harmscr(:,ii)))');
                                                        plot(EEG_FPCA_VARMAX.harmscr(:,ii),EEG_FPCA_VARMAX.harmscr(:,jj),".","MarkerSize",14,'color','black')                                           
                                                        line([0, 0], [min(EEG_FPCA_VARMAX.harmscr(:,jj))-100,max(EEG_FPCA_VARMAX.harmscr(:,jj))+100],'Color','red','LineStyle','--')
                                                        line([min(EEG_FPCA_VARMAX.harmscr(:,ii))-100,max(EEG_FPCA_VARMAX.harmscr(:,ii))+100], [0,0],'Color','red','LineStyle','--')
                                                        title(['The ', num2str(ii) ,' against ' ,num2str(jj), ' Varimax Harmonic scores function of ', SelectedChanLabels ])
                                                        xlabel(['Scores ', num2str(ii), ' of ',SelectedChanLabels, '( ',num2str(round(varExp(ii),3)*100),' %)']) 
                                                        ylabel(['Scores ', num2str(jj), ' of ',SelectedChanLabels, '( ',num2str(round(varExp(jj),3)*100),' %)'])
                                                        axis([ (min(EEG_FPCA_VARMAX.harmscr(:,ii))-100) (max(EEG_FPCA_VARMAX.harmscr(:,ii))+100) (min(EEG_FPCA_VARMAX.harmscr(:,jj))-100) (max(EEG_FPCA_VARMAX.harmscr(:,jj))+100) ])
                                                        text(EEG_FPCA_VARMAX.harmscr(:,ii),EEG_FPCA_VARMAX.harmscr(:,jj),strName,'Color', 'blue','FontSize',8) 
                                                        end 
                                                     end
                                                  end
                                              end
                          end             
                   end
            end
end    