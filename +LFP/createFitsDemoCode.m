function [fitresult, gof] = createFitsDemoCode(t, y, y2, shiftY, y3)
    %CREATEFITS(T,Y,Y2,SHIFTY,Y3)
    %  Create fits.
    %
%      t = (1:8)';
%      y = 10*log10(squeeze(MeanArray(1,5,1,:)));
%      y2 = 10*log10(squeeze(MeanArray(1,4,1,:)));
%      [~, col] = max(y)
%      shiftY = circshift(y, 4-col);
%      y3 = circshift(shiftY, -3);
%      [fitresult, gof] = Analysis.DelayedReach.LFP.createFitsDemoCode(t, y, y2, shiftY, y3)
    
    
    %  Data for 'untitled fit 1' fit:
    %      X Input : t
    %      Y Output: y
    %  Data for 'untitled fit 2' fit:
    %      X Input : t
    %      Y Output: y2
    %  Data for 'untitled fit 4' fit:
    %      X Input : t
    %      Y Output: shiftY
    %  Data for 'untitled fit 5' fit:
    %      X Input : t
    %      Y Output: y3
    %  Output:
    %      fitresult : a cell-array of fit objects representing the fits.
    %      gof : structure array with goodness-of fit info.
    %
    %  See also FIT, CFIT, SFIT.
    
    %  Auto-generated by MATLAB on 29-Jan-2018 17:11:01
    
    %% Initialization.
    
    % Initialize arrays to store fits and goodness-of-fit.
    fitresult = cell( 4, 1 );
    gof = struct( 'sse', cell( 4, 1 ), ...
        'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );
    
    %% Fit: 'untitled fit 1'.
    [xData, yData] = prepareCurveData( t, y );
    
    % Set up fittype and options.
    ft = fittype( 'gauss1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0];
    opts.Normalize = 'on';
    opts.StartPoint = [27.2232278615067 -0.612372435695795 2.38279502064522];
    
    % Fit model to data.
    [fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );
    
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult{1}, xData, yData );
    legend( h, 'y vs. t', 'untitled fit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel t
    ylabel y
    grid on
    
    %% Fit: 'untitled fit 2'.
    [xData, yData] = prepareCurveData( t, y2 );
    
    % Set up fittype and options.
    ft = fittype( 'gauss1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0];
    opts.Normalize = 'on';
    opts.StartPoint = [27.4202427156458 1.02062072615966 1.92329120146979];
    
    % Fit model to data.
    [fitresult{2}, gof(2)] = fit( xData, yData, ft, opts );
    
    % Plot fit with data.
    figure( 'Name', 'untitled fit 2' );
    h = plot( fitresult{2}, xData, yData );
    legend( h, 'y2 vs. t', 'untitled fit 2', 'Location', 'NorthEast' );
    % Label axes
    xlabel t
    ylabel y2
    grid on
    
    %% Fit: 'untitled fit 4'.
    [xData, yData] = prepareCurveData( t, shiftY );
    
    % Set up fittype and options.
    ft = fittype( 'gauss1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0];
    opts.Normalize = 'on';
    opts.StartPoint = [27.2232278615067 -0.204124145231932 2.14190058503196];
    
    % Fit model to data.
    [fitresult{3}, gof(3)] = fit( xData, yData, ft, opts );
    
    % Plot fit with data.
    figure( 'Name', 'untitled fit 4' );
    h = plot( fitresult{3}, xData, yData );
    legend( h, 'shiftY vs. t', 'untitled fit 4', 'Location', 'NorthEast' );
    % Label axes
    xlabel t
    ylabel shiftY
    grid on
    
    %% Fit: 'untitled fit 5'.
    [xData, yData] = prepareCurveData( t, y3 );
    
    % Set up fittype and options.
    ft = fittype( 'gauss1' );
    excludedPoints = excludedata( xData, yData, 'Indices', 1 );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0];
    opts.Normalize = 'on';
    opts.StartPoint = [26.6795968067268 -0.462910049886276 1.98447853804176];
    opts.Exclude = excludedPoints;
    
    % Fit model to data.
    [fitresult{4}, gof(4)] = fit( xData, yData, ft, opts );
    
    % Plot fit with data.
    figure( 'Name', 'untitled fit 5' );
    h = plot( fitresult{4}, xData, yData, excludedPoints );
    legend( h, 'y3 vs. t', 'Excluded y3 vs. t', 'untitled fit 5', 'Location', 'NorthEast' );
    % Label axes
    xlabel t
    ylabel y3
    grid on
    
    %% Fit: fit 5 without center and scale checked
    [xData, yData] = prepareCurveData( t, y3 );
    [idx,~,~,~] = util.outliers(y3);
    
    % Set up fittype and options.
    ft = fittype( 'gauss1' );
    excludedPoints = excludedata( xData, yData, 'Indices', idx );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0];
    opts.StartPoint = [27 4 2.4];
    opts.Exclude = excludedPoints;
    
    % Fit model to data.
    [fitresult{4}, gof(4)] = fit( xData, yData, ft, opts );
    
    % Plot fit with data.
    figure( 'Name', 'shifted fit' );
    h = plot( fitresult{4}, xData, yData, excludedPoints );
    legend( h, 'y vs. t', 'excluded points', 'shifted fit', 'Location', 'NorthEast' );
    % Label axes
    xlabel t
    ylabel y3
    grid on
    %%
     [xData, yData] = prepareCurveData( t, y );
    fs = fitresult{4};
    fs.b1 = fs.b1 + 2; %shift the fit curve to match data, only do this on a 
    % copy, as reassigning a coefficient clears the CI (and gives you a
    % warning)

    % Plot fit of shifted data with non-shifted data.
    % shifted excluded points back
    figure( 'Name', 'test fit' );
    h = plot( fs, xData, yData, excludedPoints );
    legend( h, 'y vs. t', 'excluded points', 'shifted shifted fit', 'Location', 'SouthEast' );
    % Label axes
    xlabel t
    ylabel y
    grid on
    
    % Now shifts the fit curve, but any data points that get swapped to the
    % other side are not matched with the curve on that side. (say the
    % right tail of the curve fit well with 2 consecutive points. If those
    % two points get shift to the left side, the curve is reflected now and
    % looks like it has large error near those points, but all points that
    % shifted well with the curve line up well). 
