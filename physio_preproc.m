%%     physio_preproc.m
% ----------------------------------------------------------------------- %
% Script for preprocessing PPG and ventilation raw data in BIDS format.   %
% Searches for accompanying MRI trigger recordings in UIM channel.        %
%                                                                         %
% Outputs:                                                                %
%                                                                         %
%   PPG:                                                                  %
%       Preprocessed scan-length PPG signal                               %
%       Peaks in PPG signal (+ history of adjustments)                    %
%       Peaks in derivative of PPG signal                                 %
%       Heart rate (HR)                                                   %
%       SQI: % of scan reliable                                           %
%       SQI: Qualitative QC label                                         %
%                                                                         %
%   Ventilation (RSP):                                                    %
%       Preprocessed scan-length RSP signal                               %
%       Peaks in RSP signal (+ history of adjustments)                    %
%       Troughs in RSP signal (+ history of adjustments)                  % 
%       Breathing rate (BR)                                               %
%       Respiratory flow (RF)                                             %
%       Respiration volume per time (RVT)                                 %
%       SQI: % of scan reliable                                           %
%       SQI: Qualitative QC label                                         %
%                                                                         %
%   Trigger locations detected                                            %
%                                                                         %
%   Parameters used in processing                                         %
%                                                                         %
% Recommended usage:                                                      %
%                                                                         %
% Adjust filtration and peak detection settings as needed, then visually  %
% QC. Regions where peaks need to be manually adjusted can be highlighted %
% using the brush tool. Another window will then open and allow the user  %
% to specify peak-specific search regions.                                %
%                                                                         %
%                                                                         %
% NOTE: do not manually close any figures.                                %
%                                                                         %
% Written by Mary Miedema, McGill University, 2025, adapted from          %
% github.com/mkassinopoulos/PRF_estimation/blob/master/Preprocess_Phys.m  %
% ----------------------------------------------------------------------- %

close all
clear all

% define path to BIDS directory
bids_dir = "D:\BIDS_working_dir\bids_MGH_FINAL\";

% define path to output directory
out_dir = "D:\estimation\physio_preproc\";

% manually select the subject ID
sub_num = "0059";

% define some MRI acquisition parameters
TR = 1.03;
n_slice = 60;
accel_f = 4;

% choose whether to overwrite previous QC outcomes for this subject
overwrite = false; 

% choose padding before/after scans in s
% NOTE: endpoints likely to be inaccurate, account for this as necessary
pad_time = 8;

i = 2;%5;

% convenient designation of sessions and scans
if i == 1
    session="02";
    task="rest";
elseif i == 2
    session="02";
    task="breathing";
elseif i == 3
    session="02";
    task="coldpressor";
elseif i == 4
    session="03";
    task="rest";
elseif i == 5
    session="03";
    task="breathing";
elseif i == 6
    session="03";
    task="visual";
end

% find the input physio files
bids_root = strcat("sub-",sub_num,"/ses-",session,"/func/sub-",sub_num,"_ses-",session,"_task-",task);
json_path = strcat(bids_dir, bids_root,"_physio.json");
physio_path = strcat(bids_dir, bids_root,"_physio.tsv.gz");

% define the output file
out_path = strcat(out_dir,bids_root);
out_physio = strcat(out_path,'_physio_and_triggers.mat');
out_QC = strcat(out_path,'_physio_QC.mat');
out_repro = strcat(out_path,'_run_repro.mat');

% check if QC and processing has already happened for this file
if isfile(out_QC) && isfile(out_physio) && overwrite == false
    disp(strcat("QC was previously completed; results in ",out_path))

% continue if not previously completed or overwriting previous    
else

    % check if the raw data exists
    if isfile(json_path) && isfile(physio_path)

        % temporarily copy, unzip, and load in the raw data
        copyfile(physio_path,'./temp_physio.tsv.gz')
        gunzip('./temp_physio.tsv.gz')
        t = readtable('./temp_physio.tsv', "FileType","text",'Delimiter', '\t');
        physio_array = table2array(t);
        delete './temp_physio.tsv'*           

        % load in the json formatted sidecar file
        fid = fopen(json_path); 
        raw = fread(fid,inf); 
        str = char(raw'); 
        fclose(fid); 
        physio_val = jsondecode(str);
    
        % extract the variables needed to process scan
        fs = physio_val.SamplingFrequency;
        % check if the sampling frequency is consistent
        fs_check = round(1/(physio_array(2,1)-physio_array(1,1)));
        if fs ~= fs_check
            disp('Warning: check sampling frequency.')
            disp('Using sampling frequency extracted from time array.')
            fs = fs_check;
        end
        cols = physio_val.Columns;

        % find time, triggers, ppg, and respiration data
        k_trig = find(contains(cols,'UIM'));
        %k_trig = find(contains(cols,'CBLCFMA')); (adjust if needed)
        k_ppg = find(contains(cols,'PPG'));
        k_rsp = find(contains(cols,'RSP'));

        time = physio_array(:,1);
        trig = zscore(physio_array(:,k_trig));
        ppg_raw = physio_array(:,k_ppg);
        rsp_raw = physio_array(:,k_rsp);

        % check the number of triggers and truncate the data 
        [trig_pks,trig_i] = findpeaks(trig,'MinPeakDistance',fs*TR*0.9,'MinPeakHeight',max(trig)*0.5);
        fprintf('    Now processing sub-%s, ses-%s, %s. ',sub_num,session,task)
        fprintf('The number of triggers found is %i.\n', size(trig_pks,1))
        trig_ind = zeros(size(physio_array,1),1);
        trig_ind(trig_i)=1;

        time0_i = find(time < 0.5/fs & time > -0.5/fs);
        time_start_i = round(time0_i - pad_time*fs);
        if time_start_i < 0
            fprintf(['The beginning of the recording conflicts with the ' ...
                'selected padding - the longest possible padding has' ...
                ' been implemented.\n'])
            time_start_i = 0;
        end

        trig_end_i = trig_i(end);
        time_end_i = round(trig_end_i + pad_time*fs);
        if time_end_i > size(physio_array,1)
            fprintf(['The end of the recording conflicts with the ' ...
                'selected padding - the longest possible padding has' ...
                ' been implemented.\n'])
            time_end_i = size(physio_array,1);
        end
        
        time = time(time_start_i:time_end_i);
        trig = trig(time_start_i:time_end_i);
        trig_i = round(trig_i-time_start_i+1); 
        trig_ind = trig_ind(time_start_i:time_end_i);
        ppg_raw = zscore(ppg_raw(time_start_i:time_end_i));
        rsp_raw = zscore(rsp_raw(time_start_i:time_end_i));
        N = time_end_i-time_start_i+1;

        % prepare output folder
        if ~exist(fileparts(out_path), 'dir')
            mkdir(fileparts(out_path))
        end

        % plot basic figure of all the raw data found
        figure('Name','Raw physiological data')
        t1 = tiledlayout(size(physio_array,2)-1,1,'TileSpacing','compact');
        for i_f1 = 2:size(physio_array,2)
            nexttile(i_f1-1)
            plot(physio_array(:,1),physio_array(:,i_f1))
            ylabel(strcat('Channel ',num2str(i_f1)))
        end
        title(t1,"Raw data read-in")

        %% do basic pre-processing on PPG data
        % bandpass filter
        ppg_f1 = 0.3; ppg_f2 = 10;  
        [filt_b, filt_a] = butter(2,[ppg_f1,ppg_f2]/(fs/2));
        ppg_filt = filtfilt(filt_b,filt_a,ppg_raw);
        PPG_filter = {[ppg_f1, ppg_f2], [filt_b, filt_a]};

        f2 = figure('Name','PPG: filtering');
        t2 = tiledlayout(2,1,'TileSpacing','compact');
        ax1 = nexttile;
        plot(time,ppg_raw)
        ylabel('Raw PPG')
        ax2 = nexttile;
        plot(time,ppg_filt)
        ylabel('Bandpass-filtered PPG')
        linkaxes([ax1 ax2],'x')
        xlim([100 103])

        % check the data for presence of time-locked scanner artefact
        prompt = "Add a notch filter at the scan slice frequency? Y/N [Y]: ";
        notch_txt = input(prompt,"s");
        if isempty(notch_txt)
            notch_txt = 'Y';
        end

        if notch_txt == 'Y'
            % notch filter the data
            notch_freq = n_slice/accel_f/TR; % scanner artefact frequency
            [b_notch, a_notch] = butter(2, [notch_freq-0.25 notch_freq+0.25]./(fs/2), 'stop');
            %[b_notch, a_notch] = butter(6, [notch_freq-0.15 notch_freq+0.15]./(fs/2), 'stop');
            ppg_bandpass = ppg_filt;
            ppg_filt = filtfilt(b_notch,a_notch,ppg_bandpass);
            PPG_filter{end+1} = [b_notch,a_notch]; 

            close(f2)
            f2 = figure('Name','PPG: filtering');
            t2 = tiledlayout(3,1,'TileSpacing','compact');
            ax1 = nexttile;
            plot(time,ppg_raw)
            ylabel('Raw PPG')
            ax2 = nexttile;
            plot(time,ppg_bandpass)
            ylabel('Bandpass-filtered PPG')
            ax3 = nexttile;
            plot(time,ppg_filt)
            ylabel('Notch and bandpass-filtered PPG')
            linkaxes([ax1 ax2 ax3],'x')
        end

        %% do basic pre-processing on respiratory data
        % de-trend
        filloutliers_window = 0.3*fs;
        filloutliers_ThresholdFactor = 0.3; 
        rsp_proc = detrend(rsp_raw,'linear');
        rsp_proc = filloutliers(rsp_proc,'linear','movmedian',filloutliers_window,'ThresholdFactor',filloutliers_ThresholdFactor);
        % bandpass filter
        if task == "breathing"
            rsp_f1 = 0.01; rsp_f2 = 4;
        else
            rsp_f1 = 0.01; rsp_f2 = 5;
        end
        [filt_b, filt_a] = butter(2,[rsp_f1,rsp_f2]/(fs/2));
        rsp_filt = filtfilt(filt_b,filt_a,rsp_proc);
        RSP_filter = {[rsp_f1, rsp_f2], [filt_b, filt_a]};

        %% downsample the data where necessary
        time_orig = time;
        DS_fs = [0,fs,nan];
        if fs > 250
            n_resamp = floor(fs/250); % get as close to 250 as possible
            time = downsample(time,n_resamp);
            ppg_raw = downsample(ppg_raw,n_resamp);
            ppg_filt = downsample(ppg_filt,n_resamp);
            rsp_filt = downsample(rsp_filt,n_resamp);
            fs = fs/n_resamp;
            trig_ind = zeros(size(time));
            for i = 1:size(trig_i,1)
                [~,index] = min(abs(time - time_orig(trig_i(i))));
                trig_ind(index) = 1;
            end
            disp(strcat("Data was downsampled to ",num2str(fs)," Hz."))
            DS_fs = [1,fs*n_resamp,fs];
        end

        time_10 = time(1):0.1:time(end);
        % find corresponding trigger indices
        trig_ind_10 = zeros(size(time_10));
        for i = 1:size(trig_i,1)
            [~,index] = min(abs(time_10 - time_orig(trig_i(i))));
            trig_ind_10(index) = 1;
        end        
        rsp_10 = interp1(time,rsp_filt,time_10);

        %% perform first-pass peak detection for PPG

        % find peaks
        min_peak_p = 0.008;% set manually as necessary, suggested 0.007-0.01
        min_peak_time = 0.65; % set a minimum distance between peaks in s, suggested 0.55-0.9
        min_peak_dist = fs*min_peak_time;
        

        % take derivative of the PPG signal
        ppg_d = diff(ppg_filt);
        [ppg_d_pks,ppg_d_i] = findpeaks(ppg_d,'MinPeakDistance',min_peak_dist,'MinPeakProminence',min_peak_p);
        ppg_d_pks_t = time(ppg_d_i);

        % find first zero-crossings after peak
        tol_0 = std(ppg_d)/3;
        ppg_d_0 = find(abs(ppg_d) < tol_0);
        ppg_pks = [];
        ppg_pks_t = [];
        for d_pk = 1:numel(ppg_d_i)
            d_0 = find(ppg_d_0>ppg_d_i(d_pk),1);
            zer_pk = ppg_d_0(d_0);
            ppg_pks_t = [ppg_pks_t,time(zer_pk)];
            ppg_pks = [ppg_pks,ppg_filt(zer_pk)];
        end
        PPG_auto_peak_params = {min_peak_p,min_peak_dist,tol_0};

        % EXCEPTIONAL: can use derivative peaks if HR otherwise unstable
        %ppg_pks = ppg_filt(ppg_d_i);
        %ppg_pks_t = ppg_d_pks_t;

        % calculate the HR
        HR = 60./diff(ppg_pks_t);
        HR_t = (ppg_pks_t(2:end)-ppg_pks_t(1:end-1))/2+ppg_pks_t(1:end-1);

        % plot an overview
        f3 = figure('Name','PPG: QC');
        t3 = tiledlayout(3,1,'TileSpacing','compact');
        ax1 = nexttile;
        p1 = plot(time,ppg_raw);
        ylabel('Raw PPG')
        ax2 = nexttile;
        p2 = plot(time,ppg_filt);
        hold on
        p3 = plot(ppg_pks_t,ppg_pks,'o');
        hold off
        ylabel('Filtered PPG')
        ax3 = nexttile;
        plot(HR_t,HR)
        title(sprintf('Heart rate (%d±%d bpm)',round(mean(HR)),round(std(HR))))
        ylabel('HR (bpm)')
        linkaxes([ax1 ax2 ax3],'x')
        xlim([time(1),time(end)])       

        %% perform manual inspection of PPG peak detection:
        
        % interactively use the brush tool to highlight
        %    1. regions of artifact which cannot be recovered (if any)
        %    2. misidentified peaks
        %    3. areas in which to add peaks manually
        % 1. should be highlighted on the top plot
        % 2. & 3. should be highlighted on the middle plot

        % Start brushing mode and wait for user to hit "Enter" when done
        brush on
        disp('Press Enter when finished QC for PPG (Figure 3).')
        pause

        check_lines = {p1, p2 ,p3};
        var_names = {'ppg_bad_data' 'ppg_pk_search' 'ppg_wrong_pks'};

        % save brushed areas to new variables
        for k = 1:numel(var_names)
            assignin('base',var_names{k},get(check_lines{k}, 'BrushData'))
        end

        % remove any wrongly identified PPG peaks
        ppg_wrong_ind = find(ppg_wrong_pks);
        ppg_pks(ppg_wrong_ind) = [];
        ppg_pks_t(ppg_wrong_ind) = [];

        %% manually add in new peaks where detection failed

        % find each brushed region from Figure 3
        search_starts = strfind([0 ppg_pk_search],[0 1]);
        search_ends = strfind([ppg_pk_search 0],[1 0]);

        % for each region
        n_stitch = 0;
        PPG_added_pks = [];
        for s = 1:size(search_starts,2)

            % extend plot for context (if possible)
            context_pad = 1;

            s1 = search_starts(s) - context_pad*fs;
            s2 = search_ends(s) + context_pad*fs;
            if s1 < 1
                s1 = 1;
            end
            if s2 > N-1
                s2 = N-1;
            end
            % if brush selected several small regions, join them
            if s < size(search_starts,2)
                if s2 > search_starts(s+1)
                    n_stitch = n_stitch + 1;
                    continue
                else
                    n_stitch = 0;
                    
                end
            end

            s1 = round(search_starts(s-n_stitch) - context_pad*fs);
            s2 = round(s2);
            s_cont = search_starts(s-n_stitch);

            % get correct peaks within this time interval
            pre_pks_i = find(ppg_pks_t > time(s1) & ppg_pks_t < time(s2));
            pre_pks_t = ppg_pks_t(pre_pks_i);
            pre_pks = ppg_pks(pre_pks_i);

            f4 = figure('Name','PPG: place new peaks');
            % zoom in on and highlight the region(s) in which the peak should be placed
            t4 = tiledlayout(3,1,'TileSpacing','compact');
            ax1 = nexttile([2 1]);
            ps = plot(time(s1:s2),ppg_filt(s1:s2));
            hold on
            plot(pre_pks_t,pre_pks,'o');
            xregion(time(s_cont),time(search_ends(s)),FaceColor="g")
            ylabel('Processed PPG')
            ax2 = nexttile;
            plot(time(s1:s2),ppg_d(s1:s2))
            ylabel('Derivative of processed PPG')
            linkaxes([ax1 ax2],'x')

            brush on
            disp(['Press Enter when finished selecting region(s) to search ' ...
                'for peak (Figure 4).'])
            pause

            assignin('base','ps_i',get(ps, 'BrushData'))

            % break into separate peak searches

            search_starts_i = strfind([0 ps_i],[0 1]);
            search_ends_i = strfind([ps_i 0],[1 0]);

            for s_i = 1:size(search_starts_i,2)

                s1_i = search_starts_i(s_i) + s1;
                s2_i = search_ends_i(s_i) + s1;

                % find peak in selected region
                [~,max_s_i] = max(ppg_d(s1_i:s2_i));
                d_0 = find(ppg_d_0>(max_s_i+s1_i),1);
                zer_pk = ppg_d_0(d_0);
                ppg_pks_t = [ppg_pks_t,time(zer_pk)];
                ppg_pks = [ppg_pks,ppg_filt(zer_pk)];

                PPG_added_pks = [PPG_added_pks,time(zer_pk)];

            end
            close(f4)
        end

        % sort the list of peaks to be in chronological order
        [ppg_pks_t,sort_pks] = sort(ppg_pks_t);
        ppg_pks = ppg_pks(sort_pks);

        % add new peaks based on linear interpolation of the heart rate
        interp_starts = strfind([0 ppg_bad_data],[0 1]);
        interp_ends = strfind([ppg_bad_data 0],[1 0]);

        PPG_interp_pks = [];

        for s = 1:size(interp_starts,2)
            t1 = time(interp_starts(s));
            t2 = time(interp_ends(s));

            % find any peaks detected in bad data
            pk1 = find(ppg_pks_t > t1,1)-1;
            pk2 = find(ppg_pks_t > t2,1);
            bad_pks_ind = find(ppg_pks_t > t1 & ppg_pks_t < t2);

            % find the nearest reliably detected peaks
            interp_n = 5;
            % TO DO: don't break things if bad section is too close to
            % start or end of recording
            mean_diff = mean(diff(ppg_pks_t(pk1-interp_n:pk1)) ...
                + diff(ppg_pks_t(pk2:pk2+interp_n)))/2;

            % remove the bad peaks
            ppg_pks(bad_pks_ind) = [];
            ppg_pks_t(bad_pks_ind) = [];
            
            % add interpolated new peaks
            n_new = round((ppg_pks_t(pk1+1) - ppg_pks_t(pk1))/mean_diff);
            %n_diff = (ppg_pks_t(pk2) - ppg_pks_t(pk1))/n_new;
            guess_pks = linspace(ppg_pks_t(pk1),ppg_pks_t(pk1+1),n_new+1);

            for s_i = 2:n_new
                % find the nearest possible point
                t_i = find(time>guess_pks(s_i),1);
                ppg_pks_t = [ppg_pks_t time(t_i)];
                ppg_pks = [ppg_pks ppg_filt(t_i)];
                PPG_interp_pks = [PPG_interp_pks time(t_i)];
            end
          
        end

        % re-sort peaks chronologically
        [ppg_pks_t,sort_pks] = sort(ppg_pks_t);
        ppg_pks = ppg_pks(sort_pks);

        % recalculate HR
        HR = 60./diff(ppg_pks_t);
        HR_t = (ppg_pks_t(2:end)-ppg_pks_t(1:end-1))/2+ppg_pks_t(1:end-1);
        HR_10 = interp1(HR_t,HR,time_10);

        % generate & save interactive figure with final peaks
        close(f3)
        fo = figure('Name','PPG: overview');
        t3 = tiledlayout(3,1,'TileSpacing','compact');
        ax1 = nexttile;
        p1 = plot(time,ppg_raw);
        ylabel('Raw PPG (V)')
        ax2 = nexttile;
        p2 = plot(time,ppg_filt);
        hold on
        p3 = plot(ppg_pks_t,ppg_pks,'o');
        hold off
        ylabel('Filtered PPG  (a.u.)')
        ax3 = nexttile;
        plot(HR_t,HR)
        title(sprintf('Heart rate (%d±%d bpm)',round(mean(HR)),round(std(HR))))
        ylabel('HR (bpm)')
        linkaxes([ax1 ax2 ax3],'x')
        xlim([time(1),time(end)])
        savefig(strcat(out_dir,bids_root,'_ppg.fig'))
        disp(strcat("QC, peak detection, and HR calculation completed for PPG data."))
        qc_prompt = "Out of the following options, \n" + ...
            "1 - Good, minimal corrections required \n" + ...
            "2 - Fine, usable with corrections \n" + ...
            "3 - Poor, not usable for most purposes \n" + ...
            "4 - Some raw data missing \n" + ...
            "5 - Otherwise problematic, consider further \n" + ...
            "Assign a qualitative label to this PPG data: ";
        PPG_q_qual = input(qc_prompt,"s");

        close(f2)
        close(fo)

        %% perform first-pass peak detection for ventilation data
        % find peaks
        min_peak_p = 0.2; % set manually as necessary
        min_peak_time = 2; % set a minimum distance between peaks in s
        min_peak_dist = fs*min_peak_time;

        [rsp_pks,rsp_i] = findpeaks(rsp_filt,'MinPeakDistance',min_peak_dist,'MinPeakProminence',min_peak_p);
        rsp_pks_t = time(rsp_i);
        [rsp_trs,rsp_tr_i] = findpeaks(-rsp_filt,'MinPeakDistance',min_peak_dist,'MinPeakProminence',min_peak_p);
        rsp_trs_t = time(rsp_tr_i);
        rsp_trs = -rsp_trs;

        RSP_auto_peak_params = {min_peak_p,min_peak_dist};

        %% correct peak/trough detection and detect bad data

        % calculate the envelopes
        respUpp = interp1([time(1),rsp_pks_t',time(end)],[rsp_pks(1),rsp_pks',rsp_pks(end)],time);
        respLow = interp1([time(1),rsp_trs_t',time(end)],[rsp_trs(1),rsp_trs',rsp_trs(end)],time);
        

        f5 = figure('Name','Ventilation: QC');
        t5 = tiledlayout(3,1,'TileSpacing','compact');
        ax1 = nexttile;
        r1 = plot(time,rsp_filt);
        ylabel('Ventilation (a.u.)')
        ax2 = nexttile([2 1]);
        r2 = plot(time,rsp_filt);
        ylabel('Ventilation (a.u.)')
        hold on
        fill([time; flipud(time)], [respLow; flipud(respUpp)], 'b', EdgeColor='none',FaceAlpha='0.15');
        r3 = plot(rsp_pks_t,rsp_pks,'o');
        r4 = plot(rsp_trs_t,rsp_trs,'o');
        hold off
        linkaxes([ax1 ax2],'x')
        xlim([time(1),time(end)])

        % Start brushing mode and wait for user to hit "Enter" when done
        brush on
        disp(['Select unusable data in the top plot. Select regions to correct' ...
            ' peaks and troughs in the bottom plot.'])
        disp('Press Enter when finished QC for ventilation (Figure 4).')
        pause

        check_lines = {r1, r2 ,r3, r4};
        var_names = {'rsp_bad_data' 'rsp_pk_search' 'rsp_wrong_pks' 'rsp_wrong_trs'};

        % save brushed areas to new variables
        for k = 1:numel(var_names)
            assignin('base',var_names{k},get(check_lines{k}, 'BrushData'))
        end

        % remove any wrongly identified peaks and troughs
        rsp_wrong_ind = find(rsp_wrong_pks);
        rsp_pks(rsp_wrong_ind) = [];
        rsp_pks_t(rsp_wrong_ind) = [];
        rsp_wrong_ind = find(rsp_wrong_trs);
        rsp_trs(rsp_wrong_ind) = [];
        rsp_trs_t(rsp_wrong_ind) = [];

        %% manually add in new peaks where detection failed

        % find each brushed region from QC
        search_starts = strfind([0 rsp_pk_search],[0 1]);
        search_ends = strfind([rsp_pk_search 0],[1 0]);

        % for each region
        n_stitch = 0;
        RSP_added_pks = [];
        RSP_added_trs = [];
        for s = 1:size(search_starts,2)

        % correct the peak detection
        % extend plot for context (if possible)
            context_pad = 3;

            s1 = search_starts(s) - context_pad*fs;
            s2 = search_ends(s) + context_pad*fs;
            if s1 < 1
                s1 = 1;
            end
            if s2 > N-1
                s2 = N-1;
            end
            % if brush selected several small regions, join them
            if s < size(search_starts,2)
                if s2 > search_starts(s+1)
                    n_stitch = n_stitch + 1;
                    continue
                else
                    n_stitch = 0;
                end
            end

            s1 = search_starts(s-n_stitch) - context_pad*fs;
            s_cont = search_starts(s-n_stitch);

            % get any correct peaks and troughs within this time interval
            pre_pks_i = find(rsp_pks_t > time(round(s1)) & rsp_pks_t < time(round(s2)));
            pre_pks_t = rsp_pks_t(pre_pks_i);
            pre_pks = rsp_pks(pre_pks_i);
            pre_trs_i = find(rsp_trs_t > time(round(s1)) & rsp_trs_t < time(round(s2)));
            pre_trs_t = rsp_trs_t(pre_pks_i);
            pre_trs = rsp_trs(pre_pks_i);

            f4 = figure('Name','Ventilation: place new peaks and troughs');
            % zoom in on and highlight the region(s) in which the peak should be placed
            ps = plot(time(s1:s2),rsp_filt(s1:s2));
            hold on
            plot(pre_pks_t,pre_pks,'o');
            plot(pre_trs_t,pre_trs,'o');
            xregion(time(s_cont),time(search_ends(s)),FaceColor="g")
            ylabel('Processed Ventilation')

            brush on
            disp(['Press Enter when finished selecting region(s) to search ' ...
                'for peaks and troughs (Figure 4).'])
            pause

            assignin('base','ps_i',get(ps, 'BrushData'))

            % break into separate peak searches

            search_starts_i = strfind([0 ps_i],[0 1]);
            search_ends_i = strfind([ps_i 0],[1 0]);

            for s_i = 1:size(search_starts_i,2)

                s1_i = search_starts_i(s_i) + s1;
                s2_i = search_ends_i(s_i) + s1;

                % look for whatever is here
                [~,max_s_i] = max(abs(rsp_filt(s1_i:s2_i)));
                % check if found a peak or trough
                max_i = round(s1_i + max_s_i);
                r_val = rsp_filt(max_i);
                if mean(diff(diff(rsp_filt(s1_i:s2_i)))) > 0
                    % found a trough
                    rsp_trs = [rsp_trs;rsp_filt(max_i)];
                    rsp_trs_t = [rsp_trs_t;time(max_i)];
                    RSP_added_trs = [RSP_added_trs;time(max_i)];
                else
                    % found a peak
                    rsp_pks = [rsp_pks;rsp_filt(max_i)];
                    rsp_pks_t = [rsp_pks_t;time(max_i)];
                    RSP_added_pks = [RSP_added_pks;time(max_i)];
                end                

            end
            close(f4)
        end

        % sort the lists of peaks and troughs to be in chronological order
        [rsp_pks_t,sort_pks] = sort(rsp_pks_t);
        rsp_pks = rsp_pks(sort_pks);
        [rsp_trs_t,sort_pks] = sort(rsp_trs_t);
        rsp_trs = rsp_trs(sort_pks);

        %% flag any problem areas for respiratory phase calculation

        % create a histogram like the one that will be used for finding the
        % phase
        NB=500;

        rsp_hist = smooth(rsp_filt, 1*fs);
        figure('Name','Ventilation: phase QC')
        t6 = tiledlayout(3,1,'TileSpacing','compact');
        nexttile
        histogram(rsp_hist,NB)
        title('Distribution of ventilation measurements')
        ylabel('Counts')
        xlabel('Ventilation amplitude')
        % calculate respiratory phase, then
        % highlight any sections of the data where respiratory phase should
        % be interpolated
        NT=length(rsp_hist);
        Phi=zeros(NT,1);        
        resp_der=diff(rsp_hist); resp_der=[resp_der;resp_der(end)];
        resp_der = filloutliers(resp_der,'linear','movmedian',0.5*fs);
        NB=1000;
        [Val,edges]=histcounts(rsp_hist,NB);
        sign_prev=1;
        for i=1:NT    
            v=rsp_hist(i);    
            [~,edge]=min(abs(v-edges));
            if (edge-v)>0 
                edge=edge-1; 
            end
            area=sum(Val(1:edge));
            sign_resp=sign(resp_der(i));
            if i > 2 && i < NT-2
                if  abs(sum(abs(resp_der(i-2:i+2)))) < 0.0005
                    sign_resp=sign_prev;
                end
            end
            sign_prev = sign_resp;
            Phi(i)=pi*area*sign_resp/NT;    
        end
        nexttile([2 1])
        plot(time,rsp_hist,'-.')
        ylabel('Respiratory phase')
        hold on
        ph1 = plot(time,Phi);

        % Start brushing mode and wait for user to hit "Enter" when done
        brush on
        disp(['Press Enter when finished highlighting regions for phase ' ...
            'interpolation (Figure 5).'])
        pause

        assignin('base',"rsp_phase_interp",get(ph1, 'BrushData'))

        %% calculate other ventilation metrics and output overview figure
        % calculate breathing rate (BR)
        BR = 60./diff(rsp_trs_t);
        time_BR = [time(1);(rsp_trs_t(2:end)-rsp_trs_t(1:end-1))/2+rsp_trs_t(1:end-1);time(end)];
        BR = interp1(time_BR,[BR(1),BR',BR(end)],time_10);

        % calculate respiratory flow (RF)
        rsp_s = smooth(rsp_10,10*1.5) ;
        RF = diff(rsp_s); RF=[0;RF(:)]; RF = RF.^2;

        % calculate RVT
        % first re-calculate upper and lower envelopes
        respUpp = interp1([time(1),rsp_pks_t',time(end)],[rsp_pks(1),rsp_pks',rsp_pks(end)],time_10);
        respLow = interp1([time(1),rsp_trs_t',time(end)],[rsp_trs(1),rsp_trs',rsp_trs(end)],time_10);
        RVT = ((respUpp-respLow).*BR)';


        % plot an overview
        figure('Name','Ventilation: overview')
        t7 = tiledlayout(4,1,'TileSpacing','compact');
        ax1 = nexttile;
        plot(time,rsp_filt)
        ylabel('Filtered ventilation (a.u.)')
        hold on
        fill([time_10, fliplr(time_10)], [respLow, fliplr(respUpp)], 'b', EdgeColor='none',FaceAlpha='0.15');
        plot(rsp_pks_t,rsp_pks,'o');
        plot(rsp_trs_t,rsp_trs,'o');
        hold off
        ax2 = nexttile;
        plot(time_10,RF)
        ylabel('Respiratory flow (a.u.)')
        ax3 = nexttile;
        plot(time_10,BR)
        ylabel('Breathing rate (rpm)')
        ax4 = nexttile;
        plot(time_10,RVT)
        ylabel('Respiration volume per time (a.u.)')
        linkaxes([ax1 ax2 ax3 ax4],'x')
        xlim([time(1),time(end)])
        savefig(strcat(out_dir,bids_root,'_resp.fig'))

        % TO DO: could add more advanced correction techniques e.g.
        % search for/highlight over regions of clipping and interpolate
        % search for belt slippage

        disp(strcat("QC, peak detection, and metric computation completed for ventilation data."))
        qc_prompt = "Out of the following options, \n" + ...
            "1 - Good, minimal corrections required \n" + ...
            "2 - Fine, usable with corrections \n" + ...
            "3 - Poor, not usable for most purposes \n" + ...
            "4 - Some raw data missing \n" + ...
            "5 - Otherwise problematic, consider further \n" + ...
            "Assign a qualitative label to this RSP data: ";
        RSP_q_qual = input(qc_prompt,"s");

        %% export final physio object
        Fs = fs;
        Fs_10 = 10;
        PPGderivlocs = ppg_d_pks_t;
        PPGlocs = ppg_pks_t;
        resp = rsp_filt;
        resp_10 = rsp_10;
        save(out_physio,'trig_ind','trig_ind_10','time','time_10','Fs','TR','PPGlocs','PPGderivlocs','HR','HR_10','Fs_10','resp','resp_10','BR','RVT','RF','rsp_phase_interp')

        %% export QC physio object
        if exist('ppg_bad_data','var') == 1
            PPG_per_bad = sum(ppg_bad_data)/length(ppg_bad_data);
        else
            PPG_per_bad = 0;
        end
        if exist('rsp_bad_data','var') == 1
            RSP_per_bad = sum(rsp_bad_data)/length(rsp_bad_data);
        else
            RSP_per_bad = 0;
        end
        save(out_QC,'PPG_per_bad','PPG_q_qual','RSP_per_bad','RSP_q_qual')

        %% export reproducibility object
        
        save(out_repro,'PPG_filter','RSP_filter','DS_fs', ...
            'PPG_auto_peak_params','ppg_bad_data','ppg_pk_search','PPG_interp_pks','PPG_added_pks', ...
            'RSP_auto_peak_params','rsp_bad_data','rsp_pk_search','RSP_added_pks','RSP_added_trs')

    else
        disp(strcat("Raw data not found! Check if ",...
            physio_path,"is expected to exist."));
    end
end