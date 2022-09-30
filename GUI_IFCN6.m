function GUI_IFCN6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Graphical user interface for spike annotation to get IFCN 6 criteria
    % and predicted probability from a LR model.
    %
    % This is the code demo prepared for our manuscript CLINPH-D-22-15528: 
    % "A quantitative approach to evaluating interictal epileptiform
    % discharges based on interpretable quantitative criteria”, under
    % review with Journal of Clincal Neurophysiology.
    %
    % Developed by Dr. Jin Jing (jjing@mgh.harvard.edu)
    % 2022-Sep-30
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %%
    warning('off','all')
    addpath('./Tools/')
    
    % Figure and axis settings
    f = figure('units','normalized','position',[0 0.0333 1.0000 0.9475]);
    set(f, 'MenuBar', 'none');
    set(f, 'ToolBar', 'none');
    set(f,'HitTest','off')
    set(f,'WindowButtonDownFcn',@clicks_Callback);
    set(f,'ButtonDownFcn',@clicks_Callback);
    set(f,'KeyPressFcn',@keys_Callback);

    Ax_model    =  subplot('position',[.05 .87 .7 .1], 'Parent',  f);
    Ax_EEG      =  subplot('position',[.05 .06 .7 .8], 'Parent',  f);
    Ax_spike    =  subplot('position',[.78 .60 .21 .37], 'Parent',  f);
    Ax_features = {subplot('position',[.78 .06 .21 .05], 'Parent',  f);
                   subplot('position',[.78 .15 .21 .05], 'Parent',  f);
                   subplot('position',[.78 .24 .21 .05], 'Parent',  f);
                   subplot('position',[.78 .33 .21 .05], 'Parent',  f);
                   subplot('position',[.78 .42 .21 .05], 'Parent',  f);
                   subplot('position',[.78 .51 .21 .05], 'Parent',  f)};

    startx = uicontrol(f,'style','pushbutton','units','normalized',  'position',[0.4690    0.58    0.050    0.040], 'string','Start','fontsize',15,'callback',@fcn_start);      
    set(Ax_EEG, 'Visible', 'off'); 
    set(Ax_spike, 'Visible', 'off'); 
    set(Ax_model, 'Visible', 'off'); 
    for j = 1:6
        set(Ax_features{j}, 'Visible', 'off');  
    end

    %% Initialization
    dataDir =  './Data/';
    tmp = load('model_IFCN6');
    features = tmp.features;
    mdl = tmp.mdl;
    feature_mean = tmp.feature_mean;
    feature_std = tmp.feature_std;
    spike = tmp.spike;
    locs =  tmp.locs;
    spike_centers = tmp.spike_centers;
    w = 10;
    ssd_thr = .42;
    nFeatures = 92;
    data = [];
    Fs = [];
    timeStampo = [];
    file = '';
    zScale = 1/100;
    ii = [];
    yp1=[];
    yp2=[];
    yb2 =[];
    pp_x =[];
    pp_y = [];
    ff = [];
    montage = 'Average';
    montages = {'Monopolar', 'Average' , 'L-Bipolar'}; 
    channel_withspace_bipolar = {'Fp1-F7' 'F7-T3' 'T3-T5' 'T5-O1' '' 'Fp2-F8' 'F8-T4' 'T4-T6' 'T6-O2' '' 'Fp1-F3' 'F3-C3' 'C3-P3' 'P3-O1' '' 'Fp2-F4' 'F4-C4' 'C4-P4' 'P4-O2' '' 'Fz-Cz'  'Cz-Pz'};
    channel_withspace_average = {'Fp1-avg' 'F3-avg' 'C3-avg' 'P3-avg' 'F7-avg' 'T3-avg' 'T5-avg' 'O1-avg' '' 'Fz-avg' 'Cz-avg' 'Pz-avg' '' 'Fp2-avg' 'F4-avg' 'C4-avg' 'P4-avg' 'F8-avg' 'T4-avg' 'T6-avg' 'O2-avg'};
    channel_withspace_monopolar = {'Fp1' 'F3' 'C3' 'P3' 'F7' 'T3' 'T5' 'O1' '' 'Fz' 'Cz' 'Pz' '' 'Fp2' 'F4' 'C4' 'P4' 'F8' 'T4' 'T6' 'O2'};
    
    uiwait(f);
    
    %% main 
    while true
        
        % parse data
        t_center = ii; 
        t_o = t_center - (w/2*Fs); 
        t_1 = t_center + (w/2*Fs)-1;
        
        seg_  = data(:, max(1, t_o) : min(size(data, 2), t_1));
        seg = NaN(size(data, 1), w*Fs);
        
        yp1_seg_  = yp1(max(1, t_o) : min(size(data, 2), t_1));yp2_seg_  = yp2(max(1, t_o) : min(size(data, 2), t_1));yb2_seg_  = yb2(max(1, t_o) : min(size(data, 2), t_1));
        yp1_seg = NaN(1, w*Fs);yp2_seg = NaN(1, w*Fs);yb2_seg = NaN(1, w*Fs);
 
        if size(seg_, 2)<size(seg, 2) 
            if t_o<1 
                seg(:, end-size(seg_, 2)+1:end) = seg_;
                yp1_seg(end-size(seg_, 2)+1:end) = yp1_seg_;yp2_seg(end-size(seg_, 2)+1:end) = yp2_seg_;yb2_seg(end-size(seg_, 2)+1:end) = yb2_seg_;
   
            elseif t_1>size(data, 2)
                seg(:, 1:size(seg_, 2)) = seg_;
                yp1_seg(1:size(seg_, 2)) = yp1_seg_;yp2_seg(1:size(seg_, 2)) = yp2_seg_;yb2_seg(1:size(seg_, 2)) = yb2_seg_;
                
            end
        else
            seg = seg_;
            yp1_seg = yp1_seg_;yp2_seg = yp2_seg_;yb2_seg = yb2_seg_;
            
        end
        eeg = seg(1:19,:);
        gap = NaN(1, size(eeg, 2));
        
        switch montage
            case 'L-Bipolar'
                seg = fcn_bipolar(eeg);
                seg_disp = [seg(1:4,:); gap; seg(5:8,:); gap; seg(9:12,:); gap; seg(13:16,:); gap; seg(17:18,:)];
                channel_withspace = channel_withspace_bipolar;
                
            case 'Average' 
                seg = eeg - repmat(mean(eeg, 1), size(eeg, 1), 1);
                seg_disp = [seg(1:8,:); gap; seg(9:11,:); gap; seg(12:19,:)];
                channel_withspace = channel_withspace_average;
                
            case 'Monopolar'
                seg =  eeg;
                seg_disp = [seg(1:8,:); gap; seg(9:11,:); gap; seg(12:19,:)];
                channel_withspace = channel_withspace_monopolar;
        end
        seg_disp(seg_disp>300) = 300;
        seg_disp(seg_disp<-300)=-300;
       
        tto = t_o; tt1 = t_1; tt = tto:tt1;
        M = size(seg_disp, 1);
        DCoff = repmat(flipud((1:M)'), 1, size(seg_disp, 2));
        timeStamps = datestr(timeStampo + seconds(round(tto/Fs):2:round(tt1/Fs)), 'hh:MM:ss');
       
        % EEG panel
        set(f,'CurrentAxes',Ax_EEG); cla(Ax_EEG)
        hold(Ax_EEG, 'on')
            for iSec = 2:round((tt1-tto+1)/Fs)
                ta = tto + Fs*(iSec-1);
                line([ta ta],  [0 M+1], 'linestyle', '--', 'color', [.5 .5 .5])
            end

            plot(Ax_EEG, tt, zScale*seg_disp+DCoff, 'k', 'linewidth', 1.0);  box on;    
            set(Ax_EEG, 'ytick',1:M,'yticklabel',fliplr(channel_withspace),'box','on', 'ylim', [0 M+1],'xlim', [tto tt1+1],'xtick',round(tt(1):2*Fs:tt(end)+1),'xticklabel',timeStamps)
            
            dt = tt1-tto+1;
            a = round(dt*4/5);
            xa1 = tto+[a a+Fs-1]; ya1 = [3 3];
            xa2 = tto+[a a]; ya2 = ya1+[0 100*zScale];
            text(xa1(1)-.7*a/10,    mean(ya2), '100\muV','Color', 'b','FontSize',12);
            text(mean(xa1), 2.5, '1 sec','Color', 'b','FontSize',12, 'horizontalalignment', 'center');        
            line(xa1,ya1, 'LineWidth', 2, 'Color','b');
            line(xa2,ya2, 'LineWidth', 2, 'Color','b');
        hold(Ax_EEG, 'off')

        % ssD model panel
        set(f,'CurrentAxes',Ax_model);cla(Ax_model)
        hold(Ax_model, 'on')
            idx_r = fcn_getRisingEdge(yb2);
            iSpikes = find(idx_r<=ii, 1, 'last');

            title([strrep(file, '.mat', ''), '  #', num2str(iSpikes) , '/' ,num2str(length(idx_r)), ' spikes'], 'interpreter', 'none') 
           
            plot(Ax_model, tt, yp1_seg, 'r', 'linewidth', 1); 
            plot(Ax_model, tt, yp2_seg, 'b--', 'linewidth', 1); 
            plot(Ax_model, tt, yb2_seg, 'b', 'linewidth', 5); 
            
            plot(Ax_model, tt, ssd_thr*ones(1, length(tt)),'--', 'color', [.5 .5 .5])
            legend('SSD', 'SSD+bgr', 'Events','threshold')
            
            set(Ax_model, 'box','on', 'ylim', [0 1], 'xlim', [tto tt1], 'xtick','','xticklabel','')
        hold(Ax_model, 'off')

        if ~isempty(pp_x) && ~isempty(pp_y)
            set(f,'CurrentAxes',Ax_EEG); 
            try
                delete(scrStr)
            catch err
                disp(err)
            end
            hold(Ax_EEG, 'on')
                plot(Ax_EEG, tt(pp_x), zScale*seg_disp(pp_y, pp_x)+DCoff(pp_y, 1), 'b.', 'markersize', 20)
            hold(Ax_EEG, 'off')
        end

        if length(pp_x)==5
            p5 = pp_x;

            tmp_t1 = p5(1)-Fs*1.5+1; tmp_t2 = p5(1)+Fs*1.5;
            tmp_seg = eeg(:, tmp_t1:tmp_t2); % c2

            switch montage
            case 'L-Bipolar'
                seg = fcn_bipolar(tmp_seg);
                ch_disp = [1:4 NaN 5:8 NaN 9:12 NaN 13:16 NaN 17:18];

            case 'Average' 
                seg = tmp_seg - repmat(nanmean(tmp_seg, 1), size(tmp_seg, 1), 1);
                ch_disp = [1:8 NaN 9:11 NaN 12:19];

            case 'Monopolar'
                seg = tmp_seg;
                ch_disp = [1:8 NaN 9:11 NaN 12:19];
            end

            pp5 = p5-tmp_t1+1;
            seg_car = tmp_seg - repmat(mean(tmp_seg, 1), size(tmp_seg, 1), 1);

            pk_v = seg_car(:, pp5(2));

            ff = fcn_getFeatures(seg(ch_disp(pp_y), :), pk_v, pp5);
            score = round(100*predict(mdl, ff))/100;

            try
                delete(scrStr)
            catch err
                disp(err)
            end
            hold(Ax_EEG, 'on')
                scrStr = text(Ax_EEG, tt(pp_x(5)), DCoff(pp_y, 1)+.5, ['\color{blue}y_p=', num2str(score)], 'color', 'r', 'fontsize', 20);
            hold(Ax_EEG, 'off')
        end

        % IFCN heatmap panel
        featureNames = {'spiky', 'asymmetry', 'duration', 'slow-wave', 'outstand', 'field'};
        for kk = 1:6
            k = 7-kk;
            set(f,'CurrentAxes',Ax_features{k});cla(Ax_features{k})
            hold(Ax_features{k}, 'on')
                colormap hot
                xm = features(kk, :);
                imagesc(Ax_features{k}, 0:8, [min(xm) max(xm)], -xm)
                plot(Ax_features{k}, 0:8, xm, '-', 'color', [0.0265 0.6137 0.8135], 'linewidth', 3)
                d = max(xm)-min(xm);
                ylim([min(xm)-d*.2 max(xm)+d*.2])
                xlim([-.5 8.5])
                axis xy
                set(Ax_features{k}, 'ytick', [],'yticklabel', [])
                if kk~=6
                    set(Ax_features{k}, 'xtick', [],'xticklabel', [])
                else
                    set(Ax_features{k}, 'xtick', 0:8, 'xticklabel', {100*(0:7)/8, '100%'})
                end
                text(Ax_features{k}, -1.5, (min(xm)+max(xm))/2, featureNames{kk})

                if ~isempty(ff)
                    ff_k = ff(kk);
                    
                    yyk = (0:8);
                    xxk = xm;

                    xq = min(max(ff_k, min(xm)), max(xm));
                    yq = interp1(xxk,yyk, xq);
                    plot(Ax_features{k},[-.5, 8.5], [xq, xq], 'k--', 'linewidth', 2); 
                    if yq>4
                        plot(Ax_features{k}, yq, xq, 'w*', 'markersize', 15, 'linewidth', 2); 
                    else
                        plot(Ax_features{k}, yq, xq, 'k*', 'markersize', 15, 'linewidth', 2); 
                    end
                end
            hold(Ax_features{k}, 'off')
        end
      
        % spike example panel
        set(f,'CurrentAxes',Ax_spike);cla(Ax_spike)
        hold(Ax_spike, 'on')
        
            fill(Ax_spike, locs(3):locs(5), spike(locs(3):locs(5)), 'r', 'Edgecolor', 'w', 'facealpha',.1)
            plot(Ax_spike, spike, 'b', 'linewidth', 3)
            
            plot(Ax_spike, locs(2), spike(locs(2))+20, 'ro', 'markersize', 20, 'markerfacecolor', 'r')
            text(Ax_spike, locs(2), spike(locs(2))+20, '1', 'color', 'w', 'fontsize', 15, 'horizontalalignment', 'center')
            text(Ax_spike, locs(2)+2, spike(locs(2))+10,   'Spiky', 'fontsize', 15)

            plot(Ax_spike, locs(2)-2, spike(locs(2))-100, 'ro', 'markersize', 20, 'markerfacecolor', 'r')
            text(Ax_spike, locs(2)-2, spike(locs(2))-100, '2', 'color', 'w', 'fontsize', 15, 'horizontalalignment', 'center')
            text(Ax_spike, locs(1)-2, spike(locs(1))-10+3,   'Asymmetry', 'fontsize', 15)

            plot(Ax_spike, locs(2)-2, spike(locs(2))-160, 'ro', 'markersize', 20, 'markerfacecolor', 'r')
            text(Ax_spike, locs(2)-2, spike(locs(2))-160, '3', 'color', 'w', 'fontsize', 15, 'horizontalalignment', 'center')
            text(Ax_spike, (locs(1)+locs(3))/2-3,spike(locs(2))-180, 'Duration','fontsize', 15)

            plot(Ax_spike, locs(4), spike(locs(2))-100, 'ro', 'markersize', 20, 'markerfacecolor', 'r')
            text(Ax_spike, locs(4), spike(locs(2))-100, '4', 'color', 'w', 'fontsize', 15, 'horizontalalignment', 'center')
            text(Ax_spike, locs(4), spike(locs(2))-120,  'Slowwave','fontsize', 15)

            plot(Ax_spike, locs(5)-10, spike(locs(2))+20, 'ro', 'markersize', 20, 'markerfacecolor', 'r')
            text(Ax_spike, locs(5)-10, spike(locs(2))+20, '5', 'color', 'w', 'fontsize', 15, 'horizontalalignment', 'center')
            text(Ax_spike, locs(5)-8, spike(locs(2))+15,  'Outstand','fontsize', 15)
            
            plot(Ax_spike, locs(5)-10, spike(locs(2))-10, 'ro', 'markersize', 20, 'markerfacecolor', 'r')
            text(Ax_spike, locs(5)-10, spike(locs(2))-10, '6', 'color', 'w', 'fontsize', 15, 'horizontalalignment', 'center')
            text(Ax_spike, locs(5)-8, spike(locs(2))-15,  'Field','fontsize', 15)

            ylim([-160 100]) 

            plot(Ax_spike, [locs(1) locs(1)], [spike(locs(2))-200 spike(locs(2))-100],  'k--');  
            plot(Ax_spike, [locs(3) locs(3)], [spike(locs(2))-200 spike(locs(2))-100],  'k--');  

            line(Ax_spike, [locs(1) locs(2) locs(3)], [spike(locs(1)) spike(locs(2)) spike(locs(3))], 'linewidth', 1, 'linestyle', '--');  

            x1 = locs(1); y1 = spike(locs(1));
            x2 = locs(2); y2 = spike(locs(2));
            a = (y1-y2)/(x1-x2);
            b = y1 - a*x1;

            y3 = -10;  x3 = (y3-b)/a;
            x4 = x3-1; y4 = a*x4+b;
            line(Ax_spike, [x3 x4 x4 x3], [y3 y4 y3 y3], 'linewidth', 1, 'linestyle', '-', 'color', 'r'); 

            dx = 0.3;dy = 2;
            line([x4+dx x4+dx x4], [y3 y3-dy y3-dy], 'linewidth', 1, 'linestyle', '-', 'color', 'r');

            text(Ax_spike, x4-1.8, (y3+y4)/2, 'S_p_r');
            text(Ax_spike, (x3+x4)/2, y3+5, '1', 'HorizontalAlignment','center');
 
            x1 = locs(3); y1 = spike(locs(3));
            x2 = locs(2); y2 = spike(locs(2));
            a = (y1-y2)/(x1-x2);
            b = y1 - a*x1;

            y3 = -10;  x3 = (y3-b)/a;
            x4 = x3+1; y4 = a*x4+b;
            line(Ax_spike, [x3 x4 x4 x3], [y3 y4 y3 y3], 'linewidth', 1, 'linestyle', '-', 'color', 'r'); 

            dx = -0.3;dy = 2;
            line(Ax_spike, [x4+dx x4+dx x4], [y3 y3-dy y3-dy], 'linewidth', 1, 'linestyle', '-', 'color', 'r');
            text(Ax_spike, x4+.3, (y3+y4)/2, 'S_p_f');
            text(Ax_spike, (x3+x4)/2, y3+5, '1', 'HorizontalAlignment','center');
            
            set(Ax_spike, 'fontsize', 12, 'xlim', [10 60])
            axis off
        hold(Ax_spike, 'off')
        
        uiwait(f);
    end

    %% Callbacks
    function fcn_start(varargin) 
        set(startx,'Visible','off','Enable', 'off');
 
        file = uigetfile(dataDir,'Select an EEG...');
        loading = uicontrol(f,'style','text','units','normalized','position',[0.469 0.580 0.052 0.05],'string','Loading...','FontSize',12);
        drawnow

        tmp = load([dataDir, file]);
        data = tmp.data;
        Fs = tmp.Fs;
        timeStampo =  tmp.startTime; 
        ii = 1;

        [B1, A1] = butter(3, .5/(Fs/2), 'high'); 
        [B2, A2] = butter(3, [50-2.5, 50+2.5]/(Fs/2), 'stop');
        data = filtfilt(B1, A1, data')'; data = filtfilt(B2, A2, data')';
        eeg_car = data(1:19, :) - mean(data(1:19, :), 1);
        
        y_amp2 = (abs(eeg_car)>=300);
        y_amp2 = sum(y_amp2, 1); y_amp2(y_amp2>0) = 1; y_amp2(y_amp2==0) =  NaN;
        
        y_amp1 = (abs(eeg_car)<=5);
        y_amp1 = sum(y_amp1, 1); y_amp1(y_amp1<=18) = NaN; y_amp1(y_amp1>18) = 1;
        
        y_amp1  = fcn_smooth(y_amp1, 8, 8, 1); y_amp1(isnan(y_amp1)) = 0; 
        y_amp2  = fcn_smooth(y_amp2 , 8, 8, 1); y_amp2(isnan(y_amp2)) = 0;
        
        yp1 = tmp.yp;   
        yp2 = yp1.*(1-tmp.Y(nFeatures, :)).*(1-y_amp1).*(1-y_amp2);   
        yb2 = NaN(1, length(yp2));
        yb2(yp2>ssd_thr) = 1;
        yb2 = fcn_smooth(yb2, 8, Fs/4, 1);  
        
        clear tmp 
        delete(loading)

        set(Ax_EEG, 'Visible', 'on');    
        set(Ax_spike, 'Visible', 'on');  
        for i = 1:6
            set(Ax_features{i}, 'Visible', 'on');  
        end
        set(Ax_model, 'Visible', 'on');  

        uiresume(f);
    end

    function dataBipolar = fcn_bipolar(data)
        % Helper function to convert EEG into L-bipolar montage
        
        dataBipolar( 1,:) = data( 1,:) - data( 5,:);  
        dataBipolar( 2,:) = data( 5,:) - data( 6,:);  
        dataBipolar( 3,:) = data( 6,:) - data( 7,:); 
        dataBipolar( 4,:) = data( 7,:) - data( 8,:);  
        dataBipolar( 5,:) = data(12,:) - data(16,:);  
        dataBipolar( 6,:) = data(16,:) - data(17,:);  
        dataBipolar( 7,:) = data(17,:) - data(18,:); 
        dataBipolar( 8,:) = data(18,:) - data(19,:);  
        dataBipolar( 9,:) = data( 1,:) - data( 2,:);  
        dataBipolar(10,:) = data( 2,:) - data( 3,:);  
        dataBipolar(11,:) = data( 3,:) - data( 4,:);  
        dataBipolar(12,:) = data( 4,:) - data( 8,:);  
        dataBipolar(13,:) = data(12,:) - data(13,:); 
        dataBipolar(14,:) = data(13,:) - data(14,:);  
        dataBipolar(15,:) = data(14,:) - data(15,:);  
        dataBipolar(16,:) = data(15,:) - data(19,:); 
        dataBipolar(17,:) = data( 9,:) - data(10,:); 
        dataBipolar(18,:) = data(10,:) - data(11,:);  
    end

    function keys_Callback(f,varargin)
        % Helper function to get keyboard input
        
        key = get(f,'CurrentKey');
        switch key
            case {'control'}
                iMontage = find(ismember(montages, montage)==1);
                iMontage = iMontage+1;
                if iMontage>3
                    iMontage = 1;
                end
                montage = montages{iMontage};            
            case 'uparrow'
                zScale =  zScale*1.5;
                
            case 'downarrow'
                zScale =  zScale/1.5;
                
            case {'leftarrow'}    
                ii = ii-round(Fs*w/2);
                pp_x = [];
                pp_y = [];
 
            case {'rightarrow'}
                ii = ii+round(Fs*w/2);
                pp_x = [];
                pp_y = [];
                   
            case 'comma'
                idx_r = fcn_getRisingEdge(yb2);
                tmp = idx_r(idx_r<ii);
                
                if isempty(tmp)
                    choice = questdlg('No more SSD spike before this!', ...
                    'Warning', ...
                    'Ok','Ok');
                    switch choice
                        case 'Ok'
                    end 
                else
                    ii = tmp(end);
                end
                pp_x = [];
                pp_y = [];
              
            case 'period'
                idx_r = fcn_getRisingEdge(yb2);
                tmp = idx_r(idx_r>ii);
                
                if isempty(tmp)
                    choice = questdlg('No more SSD spike after this!', ...
                    'Warning', ...
                    'Ok','Ok');
                    switch choice
                        case 'Ok'
                    end 
                else
                    ii = tmp(1);
                end
                pp_x = [];
                pp_y = [];
        end
        uiresume(f);
    end

    function idx_r = fcn_getRisingEdge(y)
        % Helper function to get the rising edge 
        
        yy = y;
        yy(isnan(y)) = 0;
        yy(~isnan(y)) = 1;
        
        dy = diff([0, yy]);
        idx_r = find(dy==1);
    end

    function idx_f = fcn_getFallingEdge(y)
        % Helper function to get the falling edge 
        
        yy = y;
        yy(isnan(y)) = 0;
        yy(~isnan(y)) = 1;
        
        dy = diff([yy, 0]);
        idx_f = find(dy==-1);
    end

    function y = fcn_smooth(y, gap_thr, dur_thr, true_value)
        % Helper function to smooth signal
        
        idx_r = fcn_getRisingEdge(y);
        idx_f = fcn_getFallingEdge(y);

        for k1 = 1:length(idx_r)-1
            idx_off = idx_f(k1);
            idx_on  = idx_r(k1+1);
            
            dur = idx_on - idx_off-1;
            if dur<=gap_thr  
                y(idx_off:idx_on) = true_value;
            end
        end
        
        idx_r = fcn_getRisingEdge(y);
        idx_f = fcn_getFallingEdge(y);
        for k1 = 1:length(idx_r)
            idx_on  = idx_r(k1);
            idx_off = idx_f(k1);
            
            dur = idx_off - idx_on+1;
            if dur<=dur_thr  
                y(idx_on:idx_off) = NaN;
            end
        end 
    end
    
    function xx = fcn_getFeatures(x, pk_v, p5)  
        % Helper function to get IFCN6
        
        x(isnan(x)) = 0;
 
        [s_r, s_f, aa] = fcn_getSpikienessSlope(x, p5);
        xx1 = -aa;
        xx2 = abs(s_r-s_f);
        
        [spike_dur, bg_dur_mean] = fcn_getBGdur(x,  p5);
        xx3 = bg_dur_mean/(spike_dur+eps);

        xx4 = trapz(x(p5(3):p5(5))) - trapz([p5(3), p5(5)],[x(p5(3)) x(p5(5))]);

        [P_spike, P_bg] = fcn_getPower(x, p5);
        xx5 = P_spike-P_bg;

        pk_v = (pk_v-mean(pk_v))/std(pk_v);
        xx6 = -min(mean(abs(repmat(pk_v', size(spike_centers, 1), 1) - spike_centers).^2, 2));
   
        xx = [xx1, xx2, xx3, xx4, xx5, xx6];   
        xx = (xx - feature_mean)./feature_std;
    end

    function [spike_dur, bg_dur_mean] = fcn_getBGdur(x, p5)
        % Helper function to compute IFCN feature - duration
        
        p5 = sort(p5);

        [B, A] = butter(3, 10/(Fs/2), 'low');
        x = filtfilt(B, A, x);
        x = (x-nanmean(x))/(eps+nanstd(x));
       
        zc = [];
        for i = 1:length(x)-1
            yy1 = x(i);
            yy2 = x(i+1);
            if yy1*yy2<0 || yy1==0 
                zc = [zc; i-yy1/(yy2-yy1)];
            end
        end
        zc1 = zc(zc<=p5(1)); zc2 = zc(zc>=p5(end));
 
        spike_dur = (p5(3) - p5(1));
        bg_dur_mean = mean([diff(zc1); diff(zc2)]);
    end

    function [P_spike, P_bg] = fcn_getPower(x, p5)
        % Helper function to get IFCN feature - outstand 
        
        params.movingwin= [.3, 1/Fs];  
        params.tapers   = [2,  3];
        params.fpass    = [.5, 20];
        params.Fs       = Fs;  

        x(isnan(x)) = 0;
        [spec, stimes] = mtspecgram(x',params,params.movingwin,1); 
        spec = spec';

        P = mean(pow2db(spec+eps), 1);
        P = smooth(P, 10,'sgolay')';
        P(P>25) = 25; P(P<-15) = -15;

        tt_1 = p5(1)/Fs;
        tt_2 = p5(3)/Fs;

        [~, idx1] = min(abs(stimes - tt_1));
        [~, idx2] = min(abs(stimes - tt_2));

        P_spike = mean(P(idx1:idx2));
        P_bg = mean(P([1:idx1-1, idx2+1:end]));
    end

    function [s12, s23, aa] = fcn_getSpikienessSlope(x, p5)
        % Helper function to get IFCN features - spiky and assymetry
        
        p1 = p5(2)-1; p2 = p5(2); p3 = p5(2)+1;
        s12 = abs((x(p2)-x(p1))./(p2-p1));
        s23 = abs((x(p2)-x(p3))./(p2-p3));
        aa = (pi-(atan(s12)+atan(s23)));
        
        p1 = p5(1);p2 = p5(2);p3 = p5(3);
        s12 = abs((x(p2)-x(p1))./(p2-p1));
        s23 = abs((x(p2)-x(p3))./(p2-p3));
    end

    function clicks_Callback(varargin)
        % Helper function to get in-screen clicks 
        
        click_type = get(f,'selectiontype');

        xy = get(gca,'CurrentPoint');
        kx = xy(1,1);  
        ky = xy(1,2);  
        switch click_type
            case {'normal'}  
                dd =  (tt-kx).^2;
                [~, idx_x] = min(dd);
                
                dd =  (DCoff(:, 1)-ky).^2;
                [~, idx_y] = min(dd);
                
                pp_x = unique([pp_x, idx_x]); 
                if isempty(pp_y)
                    pp_y = idx_y;
                end 
                
            case {'alt'}  
                dd =  (tt(pp_x)-kx).^2;
                [~, idx_x] = min(dd);
                pp_x(idx_x) = [];
                
                if isempty(pp_x)
                    pp_y = [];
                end    
        end
        uiresume(f);
    end

    function [spect, stimes, sfreqs] = mtspecgram(varargin)
        % Computes spectrogram of EEG data using multitaper method
        
        data1=varargin{1};
        if ~isstruct(varargin{2})
            params.pad=0;
            params.Fs=varargin{2};
            params.fpass=[0 55];
            params.err=0;
            params.trialave=0;
            params.tapers=[3 5];
            movingwin=[4 1];
        else
            params=varargin{2};
            movingwin=varargin{3};
        end

        [tapers,pad,Fs1,fpass,~,trialave,params]=getparams(params);

        data1=change_row_to_column(data1);
        [N,Ch]=size(data1);
        Nwin=round(Fs1*movingwin(1));   
        Nstep=round(movingwin(2)*Fs1);  
        nfft=max(2^(nextpow2(Nwin)+pad),Nwin);
        [sfreqs,findx]=getfgrid(Fs1,nfft,fpass);  Nf=length(sfreqs);
        params.tapers=dpsschk(tapers,Nwin,Fs1);  

        winstart=1:Nstep:N-Nwin+1;
        nw=length(winstart);
        if trialave
            S = zeros(nw,Nf);
        else
            S = zeros(nw,Nf,Ch);
        end

        for n=1:2
            indx=winstart(n):winstart(n)+Nwin-1;
            datawin=detrend(data1(indx,:));
            datawin=change_row_to_column(datawin);
            N=size(datawin,1);
            taps=dpsschk(tapers,N,Fs1);
            J=mtfftc(datawin,taps,nfft,Fs1);
            J=J(findx,:,:);
            s=squeeze(mean(conj(J).*J,2)); 
            if trialave; s=squeeze(mean(s,2));end
            S(n,:,:)=s;
        end

        for n=3:nw-1
            indx=winstart(n):winstart(n)+Nwin-1;
            datawin=detrend(data1(indx,:));
            datawin=change_row_to_column(datawin);
            J=mtfftc(datawin,taps,nfft,Fs1);
            J=J(findx,:,:);
            s=squeeze(mean(conj(J).*J,2));
            if trialave; s=squeeze(mean(s,2));end
            S(n,:,:)=s;
        end

        spect=squeeze(S);

        winmid=winstart+round(Nwin/2);
        stimes=winmid/Fs1;
    end

    function [tapers,pad,Fs,fpass,err,trialave,params]=getparams(params)
        % Helper function to convert structure params to variables 
        
        if ~isfield(params,'tapers') || isempty(params.tapers)   
             params.tapers=[3 5];
        end
        if ~isempty(params) && length(params.tapers)==3 
            TW = params.tapers(2)*params.tapers(1);
            K  = floor(2*TW - params.tapers(3));
            params.tapers = [TW  K];
        end
        if ~isfield(params,'pad') || isempty(params.pad)
            params.pad=0;
        end
        if ~isfield(params,'Fs') || isempty(params.Fs)
            params.Fs=1;
        end
        if ~isfield(params,'fpass') || isempty(params.fpass)
            params.fpass=[0 params.Fs/2];
        end
        if ~isfield(params,'err') || isempty(params.err)
            params.err=0;
        end
        if ~isfield(params,'trialave') || isempty(params.trialave)
            params.trialave=0;
        end

        tapers=params.tapers;
        pad=params.pad;
        Fs=params.Fs;
        fpass=params.fpass;
        err=params.err;
        trialave=params.trialave;
    end

    function data=change_row_to_column(data)
    % Helper routine to transform 1d arrays into column vectors
    
        dtmp=[];
        if isstruct(data)
           C=length(data);
           if C==1
              fnames=fieldnames(data);
              eval(['dtmp=data.' fnames{1} ';'])
              data=dtmp(:);
           end
        else
          [N,C]=size(data);
          if N==1 || C==1
            data=data(:);
          end
        end
    end

    function [f,findx]=getfgrid(Fs,nfft,fpass)
    % Helper function that gets the frequency grid associated with a given fft based computation
 
        if nargin < 3; error('Need all arguments'); end
        df=Fs/nfft;
        f=0:df:Fs;  
        f=f(1:nfft);
        if length(fpass)~=1
           findx=find(f>=fpass(1) & f<=fpass(end));
        else
           [~,findx]=min(abs(f-fpass));
        end
        f=f(findx);
    end

    function [tapers,eigs]=dpsschk(tapers,N,Fs)
    % Helper function to calculate tapers 
    
        if nargin < 3; error('Need all arguments'); end
        sz=size(tapers);
        if sz(1)==1 && sz(2)==2
            [tapers,eigs]=dpss(N,tapers(1),tapers(2));
            tapers = tapers*sqrt(Fs);
        elseif N~=sz(1)
            error('seems to be an error in your dpss calculation; the number of time points is different from the length of the tapers');
        end
    end

    function J=mtfftc(data,tapers,nfft,Fs)
    % Multi-taper fourier transform - continuous data
    
        if nargin < 4; error('Need all input arguments'); end
        data=change_row_to_column(data);
        [NC,C]=size(data);  
        [NK,K]=size(tapers); 
        
        if NK~=NC; error('length of tapers is incompatible with length of data'); end
        tapers=tapers(:,:,ones(1,C));  
        data=data(:,:,ones(1,K));  
        data=permute(data,[1 3 2]); 
        data_proj=data.*tapers;  
        J=fft(data_proj,nfft)/Fs;   
    end
end
