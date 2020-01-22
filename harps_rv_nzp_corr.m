function [stars,Nspec,tspan,Nrv,Nrv_corr,rvc_mean_vec,rvc_mean_vec_err,Nnoaverage,Nout_star,Nout_night,ind_same_vec,rv_std_vec,rv_med_err,nights_mean,nights_mean_err,nights_mean_std,Nrv_night,rv_std_vec_corr,rv_med_err_corr,nzp,enzp,flag_nzp] = ...
    harps_rv_nzp_corr(arm,RVtype,fname,mod_shift,jd_min,jd_max,jd_jump,Nrv_min,Nrv_min_night,sigma_outlier_star,sigma_outlier_night,init_rv_var,final_rv_var,nzp_optim,init_rv_systematics,jump_corr,smooth_nzp,smooth_span_1,smooth_span_2,...
                         sa_corr,given_drift_corr,drift_corr,erase_prop,pauseflag,display_flag,showflag,saveflag,savedata,dirname,save_fmt,unselfbias,correct_out_nzp,pvalue_outlier_nzp)
% [stars,Nspec,tspan,Nrv,Nrv_corr,rvc_mean_vec,rvc_mean_vec_err,Nnoaverage,Nout_star,Nout_night,ind_same_vec,rv_std_vec,rv_med_err,nights_mean,nights_mean_err,nights_mean_std,Nrv_night,rv_std_vec_corr,rv_med_err_corr,nzp,enzp,flag_nzp] = ...
%     harps_rv_nzp_corr(arm,RVtype,fname,mod_shift,jd_min,jd_max,jd_jump,Nrv_min,Nrv_min_night,sigma_outlier_star,sigma_outlier_night,init_rv_var,final_rv_var,nzp_optim,init_rv_systematics,jump_corr,smooth_nzp,smooth_span_1,smooth_span_2,...
%                          sa_corr,given_drift_corr,drift_corr,erase_prop,pauseflag,display_flag,showflag,saveflag,savedata,dirname,save_fmt,unselfbias,correct_out_nzp,pvalue_outlier_nzp)
% Corrects small systematic effects in HARPS RVs.
% Uses the rich data files by Trifonov et al. (2020)
%
% Input:
% arm = 'all', 'pre', or 'post' (default = 'all').
% RVtype = 'drs', 'drs_sa', 'serval', or 'serval-mlc'.
% fname = full path to the file that contains full pathes to the harps RV
%         files (default = '~\HARPS\vel_files.txt').
% mod_shift = fraction of a day to subtract from bjd so that all exposures
%             of one night will have the same floor(bjd) (default = 0.2). 
% jd_min - first night of real survey (default = 2450275 = Jul 10, 1996)
% jd_max - last night of public data (default = 2458331 = May 03, 2018)
% jd_jump - 2457174; % night of HARPS upgrade (May 29, 2015)
% Nrv_min - minimum RVs/star to be included in the analysis (default = 5).
% Nrv_min_night - minimum RVs/night to correct for its average (default = 3).
% sigma_outlier_star - Nsigma for outlier rejection per star (default = 10).
% sigma_outlier_night - Nsigma for outlier rejection per night (default = 5).
% init_rv_var - initial minimum std(RV) of an RV-loud star (m/s, default = 10).
% final_rv_var - final minimum std(RV) of an RV-loud star (m/s, default = 10).
% nzp_optim - flag to use optimum NZP corrction: all NZPs with error <
%             scatter (default = 1).
% init_rv_systematics - initial guess of the systematic RV uncertainty
%                       (added to RV uncertainties, default = 1.0 m/s).
% jump_corr - flag to correct May 2015 upgrade as jump in nzps (default = 0).
% smooth_nzp - flag to smooth NZPs and use the smothed version as NZPs (default = 1).
% smooth_span_1 - smoothing span for detecting outliers (default = 50 days).
% smooth_span_2 - smoothing span for calculating the filter model (default = 5 days)
% sa_corr = flag to correct for the secular acceleration (default: serval = 0; drs = 1).
% given_drift_corr - flag to correct for the given drift value (default = 0).
% drift_corr - flag to correct global nightly drift (default = 1).
% erase_prop - flag to erase data taken after jd_max (default: pre = 0, post = 1).
% pauseflag = flag to pause on warnings (default = 0).
% display_flag = flag to display warnings on the screen (default = 0).
% showflag - flag to make plots (default = 1).
% saveflag - flag to save plots (default = 1).
% savedata - flag to save data (default = 1).
% dirname - directory to store the results (default = '~\HARPS\results\').
% save_fmt - format of the plots (default = '.fig').
% unselfbias - flag to avoid self biasing (default = 0).
% correct_out_nzp - flag to correct for outlier NZPs (the significant ones, default = 1).
% pvalue_outlier_nzp - pvalue for outlier rejection per nzp before smoothing (default = 0.001).
%
% Output:
% stars - cell vector with the names of the stars
% Nspec - number of lines in the original .dat file
% tspan - timespan of observations (per star)
% Nrv - number of finite RVs per star, in between hd_min and jd_max
% Nrv_corr - number of finite RVs per star, corrected for NZP.
% rvc_mean_vec - vector of mean RVC per star
% rvc_mean_vec_err - vector of the error of the mean RVC per star
% Nnoaverage - number of measurements with no nightly-average correction
% Nout_star - Number of outliers per-star
% Nout_night - Number of outliers per-night
% ind_same_vec - number of nights in which the star was observed more than once
% rv_std_vec - std(RV) per star (before correction)
% rv_med_err - median RV error per star (before correction)
% nights_mean - nightly means
% nights_mean_err - nightly-mean error
% nights_mean_std - nightly-mean scatter divided by sqrt(Nrv_night)
% Nrv_night - number of good RVs per night
% rv_std_vec_corr - std(RV) per star (after correction)
% rv_med_err_corr - median RV error per star (after correction)
% nzp - NZPs before smoothing (including NaNs for no-mean nights)
% enzp - NZP errors before smoothing (including NaNs for no-mean nights)
% flag_nzp - numeric flag that explains the nature of the NZP:
%           0 - good NZP (small error and not an outlier);
%           1 - outlier NZP; 
%           2 - eNZP is larger than wrms(nzp - filter); 
%           3 - not enough RV quiet stars to calculate a NZP; 
%           4 - no NZP or filter value was calculated (inside a gap); 
%
% Created: 20180609 LT
% Last modified: 20200121 LT
%

%%
tic

%% input parameters:
if nargin<1
    arm = 'pre';
    RVtype = 'serval';
elseif nargin<3
    fname = '~\HARPS\vel_files_master_Dec19.txt';
    stars_file = '~\HARPS\stars_master_Dec19.mat';
%    warning off;
    pauseflag = 0; % flag to pause on warnings
    mod_shift = 0.2; % fraction of a day to subtract from bjd so that all exposures of one night will have the same floor(bjd)
    sigma_outlier_star = 10; % Nsigma for outlier rejection per star (default = 10).
    sigma_outlier_night = 5; % Nsigma for outlier rejection per night (default = 5).
    Nrv_min = 5; % minimum RVs/star for it to be included in the analysis (default = 5).
    Nrv_min_night = 3; % minimum RVs/night to correct for its average (default = 3).
    jd_jump = 2457162; % first night of HARPS upgrade (May 19, 2015)
    init_rv_var = 10; % initial rv std of a variable star (m/s), (default = 10).
    final_rv_var = 10; % final rv std of a variable star (m/s), (default = 10).
    init_rv_systematics = 0.8; % initial guess of the systematic rv uncertainty (added to RV uncertainties).
    smooth_nzp = 1; % flag to smooth NZPs and use the smoothed version as NZPs when no nzp was calculated.
    smooth_span_1 = 21; % smoothing span for detecting outliers (days)
    smooth_span_2 = 21; % smoothing span for calculating the filter model (days)correct_out_nzp = 1;% outlier NZPs are significant ones. The default should be correcting for them.
    correct_out_nzp = 1;% outlier NZPs are significant ones. The default should be correcting for them.
    pvalue_outlier_nzp = 0.001; % pvalue for outlier rejection per nzp before smoothing (default = 0.001).
    nzp_optim = 1; % flag to use optimum NZP corrction: all NZPs with error < scatter
    correct_nan_nzp = 1; % replace nan NZP (e.g. in gaps) with the median
    drift_corr = 1; % flag to correct for the average nightly drift
    display_flag = 0; % flag to display warnings on the screen
    showflag = 1; % flag to make plots
    saveflag = 1; % flag to save plots
    savedata = 1; % flag to save data
    unselfbias = 1; % flag to avoid self biasing (NOTE: runs for ~3 hr and makes no plots + disable remove_var option).
    dirname = ['~\HARPS\' arm '\' RVtype '\results\']; % directory to store the results
    save_fmt = '.png';
    remove_var = 0; % flag to remove specific stars to check self biasing
    variables = {'GJ253'};
end
if remove_var
    dirname = [dirname 'exclude' variables{1} '\'];
end
if strcmp(RVtype,'serval')
    sa_corr = 0; % flag to correct for the secular acceleration
    given_drift_corr = 0; % flag to correct for the given drift value (in e2ds fits headers).
else
    sa_corr = 1; % flag to correct for the secular acceleration
    given_drift_corr = 0; % flag to correct for the given drift value (in e2ds fits headers).
end
if strcmp(arm,'pre')
    jd_min = 2452936; % first night of real survey (Oct 23, 2003)
    jd_max = 2457162; % last night of HARPS pre (May 19, 2015)
    jump_corr = 0; % flag to correct May 2015 upgrade as jump in nzps.
elseif strcmp(arm,'post')
    jd_min = 2457174; % first night of HARPS post (May 31, 2015)
    jd_max = 2458252; % last night of the public survey data (May 13, 2018)
    jump_corr = 0; % flag to correct May 2015 upgrade as jump in nzps.
else 
    jd_min = 2452936; % first night of real survey (Oct 23, 2003)
    jd_jump = 2457168; % middle of HARPS upgrade (May 25, 2015)
    jd_max = 2458252; % last night of the public survey data (May 13, 2018)
    jump_corr = 0; % flag to correct May 2015 upgrade as jump in nzps.
end

%% ===========================================================================
% STARTING POINT OF CALCULATING THE SPARSE MATRIX OF RV-NIGHT
%===========================================================================

% load the list of .dat files and initialize some output vectors:
targets = loadtxt(fname);
target_num = length(targets);
load(stars_file); 
rv_std_vec = nan(length(targets),1);
Nrv = rv_std_vec;
Nspec = Nrv;
Nrv_corr = Nrv;
rvc_mean_vec = Nrv;
rvc_mean_vec_err = Nrv;
rv_med_err = Nrv;
ind_same_vec = Nrv;
Nnoaverage = zeros(length(targets),1);
Nout_star = Nrv;
tspan = Nrv;

% Initialize the matrices:
a = clock;
day_str = [num2str(a(1),'%04d') num2str(a(2),'%02d') num2str(a(3),'%02d')];
% jd_now = utc2jd(a(1),a(2),a(3),a(4),a(5),a(6));
m = jd_max-jd_min+1;
rv_mat = nan(target_num,m);
err_mat = rv_mat;
jd_mat = rv_mat;

% open a .orb file:
if savedata
    orbname = [dirname 'orb_files\harps_' RVtype '_' arm '.orb'];
    fid = fopen(orbname,'w');
end

fid_empty = fopen(['~\HARPS\' RVtype '-' arm '_empty_files.txt'],'w');
fid_nan = fopen(['~\HARPS\' RVtype '-' arm '_nan_RV_files.txt'],'w');

%% collect the times, RVs, and errors into the matrices, target by target:
for t=1:target_num
    
    % load the .dat file:
    filename = targets{t};
    starname = stars{t};
    try
        M = load(filename);
    catch
        warning([filename ' does not exist. Press a key to continue']);
        if pauseflag; pause; end
        continue
    end
    
    % Empty files can be skipped:
    if isempty(M)
        if display_flag
            display([filename ' is an empty file. Press a key to continue']);
        end
        fprintf(fid_empty,'%s\n',filename);
        if pauseflag; pause; end
        continue
    end
    
    % sort the matrix and count total number of spectra:
    M = sortrows(M,1);
    bjd_tmp = M(:,1)-mod_shift;
    Nspec(t) = length(bjd_tmp);
    
    try
        M(bjd_tmp<jd_min | bjd_tmp>jd_max,:)=[]; % Remove measurements outside [jd_min,jd_max].
    catch
        if display_flag
            display([starname ' has no RVs to remove.']);
        end
    end
    
    % Empty files should be skipped:
    if isempty(M)
        if display_flag
            display([filename ' has no RVs in range.']);
        end
        if pauseflag; pause; end
        continue
    end
    
    % read the data into vectors:
    bjd = M(:,1)-mod_shift; % Subtract mod_shift day from time stamps to have all
                      % the RVs of each night with the same floor(bjd).
    if strcmp(RVtype,'serval')
        rvc = M(:,6);
        rvc_err = M(:,7);
    elseif strcmp(RVtype,'drs') % drs RVs:
        rvc = M(:,8);
        rvc_err = M(:,9);
    elseif strcmp(RVtype,'crx') % crx
        rvc = M(:,12);
        rvc_err = M(:,13);
    else
        error('Wrong RVtype');
    end
    flag_RV = M(:,22); % numeric flag
    drift = M(:,31);
    drift(isnan(drift)) = 0;
    sa = M(:,33);
    sa(isnan(sa)) = 0;
    if given_drift_corr
       rvc = rvc - drift; 
    end
    if sa_corr
       rvc = rvc - sa; 
    end
    rv_med_err(t) = nanmedian(rvc_err);

    % remove nan values:
    ind_nan = find(isnan(rvc) | isnan(rvc_err));
    if ~isempty(ind_nan)
        if display_flag
            display([filename ' contains ' num2str(sum(isnan(rvc))) '/' num2str(length(rvc)) ' NaN RVs. and ' num2str(sum(isnan(rvc_err))) '/' num2str(length(rvc_err)) ' NaN RV errors. Removed.']);
        end
        fprintf(fid_nan,'%s:\n %d/%d NaN values\n',filename,length(ind_nan),length(rvc));
        bjd(ind_nan) = [];
        rvc(ind_nan) = [];
        rvc_err(ind_nan) = [];
        flag_RV(ind_nan) = [];
        if pauseflag; pause; end
        if isempty(rvc)
            if display_flag
                display([filename ' has no good RVs in range.']);
            end
        if pauseflag; pause; end
            continue
        end
    end
    
    % Find outliers (do not rely on RV errors - they do not account for systematics yet),
    % and flagged RVs (S/N<10 or S/N>400):
    Nrv(t) = sum(isfinite(rvc));
    tspan(t) = max(bjd) - min(bjd);
    weights = rvc_err.^-2;
    rv_mean = nanwmean(rvc,weights);
    bias = 1-1/4/Nrv(t);
    rv_std = wstd(rvc,weights)/bias;
    bad_rv = find((abs(rvc-rv_mean)>sigma_outlier_star*rv_std & abs(rvc-rv_mean)>init_rv_var)...
                  | flag_RV>0);
    rv_std_vec(t) = rv_std;
    Nout_star(t) = length(bad_rv);
    
    % Write an entry to .orb file (original RVs including the outliers but excluding RVs outside [jd_min,jd_max]):
    if savedata
        fprintf(fid,'STAR: %s\n',starname);
        for i=1:length(bjd)
            fprintf(fid,'A %f %f %f\n',bjd(i)+mod_shift,0.001*rvc(i),0.001*rvc_err(i));
        end
        fprintf(fid,'END\n');
    end
    
    % Remove stellar outliers (only for the purpose of calculating the NZPs):
    bjd(bad_rv) = [];
    rvc(bad_rv) = [];
    rvc_err(bad_rv) = [];
    
    % For nights with multiple exposures, take the weighted mean RVC in that
    % night, with its scatter as an error
    nights = floor(bjd);
    columns = nights-jd_min+1;
    ind_same_tmp = 0;
    for i=1:m
        ind_same = find(abs(columns-i)==0);
        if ~isempty(ind_same) && length(ind_same)>1
            [~,ind_min_tmp] = min(rvc_err(ind_same));
            % calculate a weihted mean values:
            bjd_ind_same = nanwmean(bjd(ind_same),rvc_err(ind_same).^-2);
            rvc_ind_same = nanwmean(rvc(ind_same),rvc_err(ind_same).^-2);
            rvc_err_ind_same = max(rvc_err(ind_same(ind_min_tmp)),wstd(rvc(ind_same),rvc_err(ind_same).^-2));
            % save the corrected values instead of the minimum err value:
            bjd(ind_same(ind_min_tmp)) = bjd_ind_same;
            rvc(ind_same(ind_min_tmp)) = rvc_ind_same;
            rvc_err(ind_same(ind_min_tmp)) = rvc_err_ind_same;
            % erase the index from the values to be erased:
            ind_same(ind_min_tmp) = [];
            % erase the the rest of the values:
            columns(ind_same) = [];
            nights(ind_same) = [];
            bjd(ind_same) = [];
            rvc(ind_same) = [];
            rvc_err(ind_same) = [];
            % count the number of binned data
            ind_same_tmp = ind_same_tmp+1;
        end
    end
    ind_same_vec(t) = ind_same_tmp;
    
    % Re-calculate the weighted mean (zero-point of the star) and error:
    if sum(isfinite(rvc))
        bias = 1-1/4/Nrv(t);
        rvc_err = sqrt(rvc_err.^2+init_rv_systematics^2);
        weights = rvc_err.^-2;
        rvc_mean = nanwmean(rvc,weights);
        rvc_mean_vec(t) = rvc_mean;
        rvc_mean_err = 1/sqrt(sum(weights));
        rvc_mean_std = wstd(rvc,weights)/sqrt(Nrv(t))/bias;
        rvc_mean_err = max(rvc_mean_err,rvc_mean_std);
        rvc_mean_vec_err(t) = rvc_mean_err;

        % Write into the matrices:
        jd_mat(t,columns) = bjd;
        rv_mat(t,columns) = rvc-rvc_mean; % stellar zero-point subtracted RVs
        err_mat(t,columns) = sqrt(rvc_err.^2+rvc_mean_err^2); % the error on the stellar zero-point is co-added to the RV error,
                                                              % so that variables and/or faint stars get lower weight in NZP calculation. 
    else
        warning([starname ' has no good RVs.']);
    end
    
end

fclose(fid_empty);
fclose(fid_nan);

% save the matrices
save([dirname 'mat_files\jd_mat'],'jd_mat','rv_mat','err_mat');

% close the .orb and all_rvs files
if savedata
    fclose(fid);
end

%% =========================================================================
% STARTING POINT OF CORRECTING THE RADIAL VELOCITIES WITH SELF BIASING
%===========================================================================
flag_nzp = zeros(1,m);

if ~unselfbias
% Find variable stars:
    novar_rv = find(rv_std_vec<init_rv_var);
    rv_std_median = nanmedian(rv_std_vec(novar_rv));
    if showflag
        disp('Before correction:');
        fprintf(1,'Median std(RV) of RV-quiet stars: %4.2f\n',rv_std_median);
        fprintf(1,'std(RV) threshold for an RV_loud star: %4.2f\n',init_rv_var);
    end
    var_rv = find(rv_std_vec>init_rv_var);
    
    % Remove also hidden-variable stars:
    for t=1:target_num
        starname = stars{t}; 
        for v = 1:length(variables)
            if strcmp(variables{v},starname) && remove_var
               var_rv = [var_rv;t];
            end
        end
    end
    
    % remove multiple entries:
    var_rv = sort(var_rv);
    var_rv(diff(var_rv)==0)=[];
    
    
%% plot std(RV) histogram
if showflag
     figure(16);
     set(16,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
     hold off
     g = logspace(0,1.05*log10(max(rv_std_vec)),ceil(10*(1.05*log10(max(rv_std_vec)))));
     hb = hist(rv_std_vec(Nrv>=Nrv_min),g);
     stairs(g,hb,'r','LineWidth',2);
     hold on;
     xlim([1 10^(1.05*log10(max(rv_std_vec)))]);
     ylim([0 1.05*max(hb)]);
     set(gca,'Xscale','log');
     ylabel('Number of stars');
     xlabel('std(RV) [m/s]');
     title([date ': Histograms of std(RV) before and after correcting for the NZPs']);
     set(gca,'FontSize',16);
     grid on
end

%% Average RVs from the same night using all stars observed in that night:

% initialize some vectors:
nights_mean = nan(1,m);
nights_median = nights_mean;
nights_mean_err = nights_mean;
nights_mean_std = nights_mean;
Nrv_night = nan(1,m);
Nout_night = nan(1,m);

% load the results of stage 1 (the sparse RV-night matrices):
load([dirname 'mat_files\jd_mat.mat'],'jd_mat','rv_mat','err_mat');

% Remove the RVs of targets with less than Nrv_min exposures (only for the purpose of calculating the NZPs):
jd_mat(Nrv<Nrv_min,:) = nan;
rv_mat(Nrv<Nrv_min,:) = nan;
err_mat(Nrv<Nrv_min,:) = nan;

% Remove the RVs of variable targets (only for the purpose of calculating the NZPs):
jd_mat(var_rv,:) = nan;
rv_mat(var_rv,:) = nan;
err_mat(var_rv,:) = nan;  

% start the plot of all RVs of RV-quiet stars and the NZPs:
if showflag
    figure(10);
    set(10,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    hold off
end
for t=1:m
    night_jds = jd_mat(:,t);
    night_rvs = rv_mat(:,t);
    night_err = err_mat(:,t);
    
    % plot nightly RVs
    if showflag
        figure(10);
        hold on
        plot(night_jds-jd_min+mod_shift,night_rvs,'.c');
    end
    
    % remove outliers per night (but plot them) only if there are enough RVs:
    good_rvs_tmp = find(isfinite(night_rvs));
    Nrv_night_tmp = length(good_rvs_tmp);
    if Nrv_night_tmp>Nrv_min_night
        night_weights = night_err.^-2;
        rv_mean_night = nanwmean(night_rvs,night_weights);
        bias = 1-1/4/Nrv_night_tmp;
        rv_std_night = wstd(night_rvs,night_weights)/bias;
        bad_rv_night = find(abs(night_rvs-rv_mean_night)>sigma_outlier_night*rv_std_night & abs(night_rvs-rv_mean_night)>init_rv_var);
        if showflag
            plot(night_jds(bad_rv_night)-jd_min+mod_shift,night_rvs(bad_rv_night),'or')
        end
        night_jds(bad_rv_night) = nan;
        night_rvs(bad_rv_night) = nan;
        night_err(bad_rv_night) = nan;
        Nout_night(t) = length(bad_rv_night);
    end
    
    %==========================================================
    % Calculate the nightly zero-point RV using only good RVs:
    good_rvs = find(isfinite(night_rvs));
    Nrv_night(t) = length(good_rvs);
    bias = 1-1/4/Nrv_night(t);
    weights = night_err.^-2;
    if sum(isfinite(weights))
        nights_mean(t) = nanwmean(night_rvs(good_rvs),weights(good_rvs));
        nights_mean_err(t) = 1/sqrt(sum(weights(good_rvs)))/bias;
        nights_median(t) = wmedian(night_rvs(good_rvs),weights(good_rvs));
        nights_mean_std(t) = wstd(night_rvs(good_rvs),weights(good_rvs))/sqrt(Nrv_night(t))/bias;
        % Take the max of nights_mean_err and nights_mean_std:
        nights_mean_err(t) = max(nights_mean_err(t),nights_mean_std(t));
    end
    %==========================================================
    
end

%% display some statistics:
    weights_nzp = nights_mean_err.^-2;
    wmean_nights_mean = nanwmean(nights_mean(Nrv_night>=Nrv_min_night),weights_nzp(Nrv_night>=Nrv_min_night));
    std_nights_mean = wstd(nights_mean(Nrv_night>=Nrv_min_night),weights_nzp(Nrv_night>=Nrv_min_night));
    median_nights_mean_err = nanmedian(nights_mean_err(Nrv_night>=Nrv_min_night));
    fprintf(1,'std(NZP) = %4.2f m/s; med(eNZP) = %4.2f m/s\n',std_nights_mean,median_nights_mean_err);
    
%% find nights where NZP correction could not be applied and mark them:
bad_mean = find(Nrv_night<Nrv_min_night);
nights_vec = 1:m;

if showflag
    figure(10)
    h = errorbar(nights_vec(isfinite(nights_mean))-0.3,nights_mean(isfinite(nights_mean)),nights_mean_err(isfinite(nights_mean)),'.k','MarkerSize',15,'LineWidth',1);
% Make cap-lines mult times long
    mult = 0;                               % elongation factor
    b = h.Bar;                              % hidden property of h=errorbar(X,Y,E)
    drawnow                                 % populate b's properties
    vd = b.VertexData;
    X = nights_vec(isfinite(nights_mean));
    N = numel(X);                           % number of error bars
    capLength = vd(1,2*N+2,1) - vd(1,1,1);  % assumes equal length on all
    newLength = capLength * mult;
    leftInds = N*2+1:2:N*6;
    rightInds = N*2+2:2:N*6;
    vd(1,leftInds,1) = [X-newLength, X-newLength];
    vd(1,rightInds,1) = [X+newLength, X+newLength];
    b.VertexData = vd;
% plot the bad-mean NZPs:
    plot(nights_vec(bad_mean)-0.3,nights_mean(bad_mean),'sr');
end

%% ==========================================================
% Optimize the NZP model:

% Put NaN for no-mean nights for the smoothing
nights_mean(bad_mean) = NaN;
nights_mean_err(bad_mean) = NaN;

    nzp = nights_mean(:);
    enzp = nights_mean_err(:);
    weight = 1./enzp.^2;
    N_sigma_0 = abs(nzp./enzp); % LT 20190623: check also distance from 0.
    prob_n_0 = 2*(1-cdf('t',N_sigma_0',Nrv_night-1,0,1));
    [nzp_filt,nzp_filt_err] = wmeanfilt02(nights_vec,nzp,weight,smooth_span_1);

if smooth_nzp
    if ~jump_corr
        if correct_out_nzp
            N_sigma = abs((nzp-nzp_filt)./enzp);
            prob_n = 2*(1-cdf('t',N_sigma',Nrv_night-1,0,1));
            outlier_nzp = find(prob_n<pvalue_outlier_nzp & prob_n_0<pvalue_outlier_nzp);
            Nout = length(outlier_nzp);
            nzp_out = nzp;
            enzp_out = enzp;
            nzp(outlier_nzp) = nan;
            enzp(outlier_nzp) = nan;
            weight = 1./enzp.^2;
            [nzp_filt,nzp_filt_err] = wmeanfilt02(nights_vec,nzp,weight,smooth_span_2);
        end
        % correct all RVs with the smoothing filter
        nights_mean = nzp_filt';
        nights_mean_err = nzp_filt_err';
        % but keep outlier NZPs
        if correct_out_nzp
            nights_mean(outlier_nzp) = nzp_out(outlier_nzp);
            nights_mean_err(outlier_nzp) = enzp_out(outlier_nzp);
            nzp(outlier_nzp) = nzp_out(outlier_nzp);
            enzp(outlier_nzp) = enzp_out(outlier_nzp);
        end
        
    elseif jump_corr
        [nzp_filt_b,nzp_filt_err_b] = wmeanfilt02(nights_vec(1:(jd_jump-jd_min)),nzp(1:(jd_jump-jd_min)),weight(1:(jd_jump-jd_min)),smooth_span_1);
        [nzp_filt_a,nzp_filt_err_a] = wmeanfilt02(nights_vec(jd_jump-jd_min+1:end),nzp(jd_jump-jd_min+1:end),weight(jd_jump-jd_min+1:end),smooth_span_1);
        nzp_filt = [nzp_filt_b;nzp_filt_a];
        nzp_filt_err = [nzp_filt_err_b;nzp_filt_err_a];
        if correct_out_nzp
            N_sigma = abs((nzp-nzp_filt)./enzp);
            prob_n = 2*(1-cdf('t',N_sigma',Nrv_night-1,0,1));
            outlier_nzp = find(prob_n<pvalue_outlier_nzp);
            Nout = length(outlier_nzp);
            nzp_out = nzp;
            enzp_out = enzp;
            nzp(outlier_nzp) = nan;
            enzp(outlier_nzp) = nan;
            weight = 1./enzp.^2;
            [nzp_filt_b,nzp_filt_err_b] = wmeanfilt02(nights_vec(1:(jd_jump-jd_min)),nzp(1:(jd_jump-jd_min)),weight(1:(jd_jump-jd_min)),smooth_span_2);
            [nzp_filt_a,nzp_filt_err_a] = wmeanfilt02(nights_vec(jd_jump-jd_min+1:end),nzp(jd_jump-jd_min+1:end),weight(jd_jump-jd_min+1:end),smooth_span_2);
            nzp_filt = [nzp_filt_b;nzp_filt_a];
            nzp_filt_err = [nzp_filt_err_b;nzp_filt_err_a];
        end
        % correct all RVs with the smoothing filter
        nights_mean = nzp_filt';
        nights_mean_err = nzp_filt_err';
        % but keep outlier NZPs
        if correct_out_nzp
            nights_mean(outlier_nzp) = nzp_out(outlier_nzp);
            nights_mean_err(outlier_nzp) = enzp_out(outlier_nzp);
            nzp(outlier_nzp) = nzp_out(outlier_nzp);
            enzp(outlier_nzp) = enzp_out(outlier_nzp);
        end
        
    end
    
    % optimum NZP corrction: use the individual NZPs with error < scatter(res):
    res_rms = wstd(nzp' - nights_mean,enzp.^-2);
    ind_optim = find(enzp<res_rms);
    ind_not_optim = find(enzp>res_rms);
    if nzp_optim
        fprintf(1,'\ncorrecting for %d/%d NZPs with errors<%4.2f m/s\n',length(ind_optim),length(ind_optim)+length(ind_not_optim),res_rms);
        nights_mean(ind_optim) = nzp(ind_optim);
        nights_mean_err(ind_optim) = enzp(ind_optim);
        nights_mean_err(ind_not_optim) = res_rms;
        flag_nzp(ind_not_optim) = 2;
    end
   
    % Keep significant NZPs and flag them:
        if correct_out_nzp
            fprintf(1,'\ncorrecting for %d significant NZPs\n',Nout);
            nights_mean(outlier_nzp) = nzp_out(outlier_nzp);
            nights_mean_err(outlier_nzp) = enzp_out(outlier_nzp);
            flag_nzp(outlier_nzp) = 1;
        end
        
    % put the NZP filter as estimate and scatter as error for no-mean nights: 
        nzp_filt_est = nzp_filt(:);
        nights_mean(bad_mean) = nzp_filt_est(bad_mean);
        fprintf(1,'\ncorrecting for %d NZPs from nights with Nrv<%d by the wmean filter +/- %4.2f\n',length(bad_mean),Nrv_min_night,res_rms);
        nights_mean_err(bad_mean) = res_rms;
        flag_nzp(bad_mean) = 3;
    
    % replace 0 and NaN NZPs (where there is no information) with the mean and rms:
        ind_nan = find(nights_mean==0 | ~isfinite(nights_mean) | ~isfinite(nights_mean_err));
        Nnan = length(ind_nan);
        if correct_nan_nzp
            fprintf(1,'\ncorrecting for %d NaN (or 0) NZPs by %4.2f +/- %4.2f\n\n',Nnan,wmean_nights_mean,std_nights_mean);
            nights_mean(ind_nan) = wmean_nights_mean;
            nights_mean_err(ind_nan) = std_nights_mean;
            flag_nzp(ind_nan) = 4;
        end
    
else % Avoid smoothing the NZPs:    
        nights_mean(:) = nzp(:);
        nights_mean_err(:) = enzp(:);
end    

    % NZP plot:
        if showflag
             if nzp_optim
                 plot(nights_vec(ind_not_optim)-0.3,nzp(ind_not_optim),'or');
             end
             if correct_out_nzp
                 plot(nights_vec(outlier_nzp)-0.3,nights_mean(outlier_nzp),'om');
             end
            plot(nights_vec-0.3,nights_mean,'.g','MarkerSize',10);
        end
% ==========================================================
    
%% some design:
if showflag
    figure(10);
    xlim([-100 jd_max-jd_min+100]);
    hold on;
    plot([-100 jd_max-jd_min+100],[final_rv_var final_rv_var],'-k');
    ylim([-final_rv_var final_rv_var]);
    plot([jd_max-jd_min+100 jd_max-jd_min+100],[-final_rv_var final_rv_var],'-k');
    xlabel (['JD - ' num2str(jd_min)],'fontsize',14);
    set(10,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    ylabel('RV [m/s]');
    set(gca,'FontSize',12);
    grid on
    if saveflag
        filename = ['harps_' arm '_' RVtype '_night_means_selfbiased_' day_str];
        saveas(gcf,[dirname 'figures\' filename save_fmt]);
    end
end

%% calculate the global nightly drift
rv_corr = rv_mat;
for j=1:m;
    rv_corr(:,j)=rv_corr(:,j)-nights_mean(j);
end
jd = jd_mat(:);
rva = rv_corr(:);
err = err_mat(:);
inan = isnan(jd);
jd(inan)=[];
rva(inan)=[];
err(inan)=[];
clear inan
modbjd = mod(jd,1)-0.5; % jd is bjd-mod_shift; 
% =====================================================================
if drift_corr
   [a1,b1,da1,db1,~,~,~,~,~,Ftest1] = linfit(modbjd,rva,err,0);
   if showflag
    figure(11)
        plot(modbjd,rva,'.')
        set(11,'units','normalized','outerposition',[0.05 0.05 0.6 0.9]);
        hold on
        grid on
        plot([min(modbjd)-0.01 max(modbjd)+0.01],b1+[min(modbjd)-0.01 max(modbjd)+0.01]*a1,'-k','LineWidth',2.5)
        ylabel('RV (after NZP correction) [m/s]');
        ylim([-final_rv_var final_rv_var])
        xlim([-0.26 0.26])
        set(gca,'xTick',[-3/12 -2/12 -1/12 0 1/12 2/12 3/12]);
        set(gca,'xTickLabel',[-6 -4 -2 0 2 4 6]);
        xlabel('t_{mid} [hr]');
        legend(['p(F_{test})_{value}  = ' num2str(Ftest1,'%3.1s') ' (' num2str(length(rva)) ' points)'],...
               ['line slope: a = ' num2str(a1./24,'%5.3f') '+/-' num2str(da1./24,'%5.3f') ' m/s/h']);
        set(gca,'FontSize',16)
        if saveflag
           filename = ['harps_' arm '_' RVtype '_nightly_drift_selfbiased_' day_str];
           saveas(gcf,[dirname 'figures\' filename save_fmt]);
        end
   end
end
% =====================================================================


%% apply the nightly correction, add its uncertainties, and write a new .orb, and avc.dat files:
if savedata
    orbname = [dirname 'orb_files\harps_' arm '_' RVtype '_avc_selfbiased_' day_str '.orb'];
    fid = fopen(orbname,'w');
    allname = [dirname 'txt_files\harps_' arm '_' RVtype '_all_rvs_selfbiased_' day_str '.tsv'];
    fid_all = fopen(allname,'w');
end
rv_std_vec_corr = nan(length(targets),1);
rv_med_err_corr = nan(length(targets),1);

for t=1:target_num
   
    % Re-load all the RVs from the dat file for correction:
    filename = targets{t};
    starname = stars{t};
    try
        M = load(filename);
    catch
        warning([filename ' does not exist. Press a key to continue']);
        if pauseflag; pause; end
        continue
    end
    
    % Empty files can be skipped:
    if isempty(M)
        if display_flag
            display([filename ' is an empty file. Press a key to continue']);
        end
        if pauseflag; pause; end
        continue
    end
    
    % sort the matrix:
    M = sortrows(M,1);
    bjd_tmp = M(:,1)-mod_shift;
    
    % read the data into vectors:
    bjd = M(:,1)-mod_shift; % Subtract mod_shift day from time stamps to have all
                      % the RVs of each night with the same floor(bjd).
    if strcmp(RVtype,'serval')
        rvc = M(:,6);
        rvc_err = M(:,7);
    elseif strcmp(RVtype,'drs') % drs RVs:
        rvc = M(:,8);
        rvc_err = M(:,9);
    elseif strcmp(RVtype,'crx') % crx
        rvc = M(:,12);
        rvc_err = M(:,13);
    end
    drift = M(:,31);
    drift(isnan(drift)) = 0;
    sa = M(:,33);
    sa(isnan(sa)) = 0;
    if given_drift_corr
       rvc = rvc - drift; 
    end
    if sa_corr
       rvc = rvc - sa; 
    end
    
    inrange = find(bjd_tmp>jd_min & bjd_tmp<jd_max);
    outrange = find(bjd_tmp<jd_min | bjd_tmp>jd_max);
    
    % Empty files should be skipped:
    if isempty(inrange)
        if display_flag
            display([filename ' has no RVs in range.']);
        end
        if pauseflag; pause; end
        continue
    end
    
    % remove values out of range:
         bjd(outrange) = [];
         rvc(outrange) = [];
         rvc_err(outrange) = [];
    
    modbjd_t = mod(bjd,1)-0.5;
    corr_t = nan(size(bjd));
    corr_err = nan(size(bjd));
    
    %==========================================================
    % subtract the nightly mean and co-add its error:
    nights = floor(bjd);
    columns = nights-jd_min+1;
    for i=1:length(bjd)
        if drift_corr
            corr_t(i) = nights_mean(columns(i)) + (b1+modbjd_t(i)*a1);
        else
            corr_t(i) = nights_mean(columns(i));
        end
        rvc(i) = rvc(i) - corr_t(i);
        if drift_corr
            corr_err(i) = sqrt(nights_mean_err(columns(i))^2 + db1^2 + (modbjd_t(i)*da1)^2);
        else
            corr_err(i) = nights_mean_err(columns(i));
        end
        rvc_err(i) = sqrt(rvc_err(i)^2 + corr_err(i)^2);
    end
    %==========================================================
     
    % for serval RVs of pre or post separately, remove again the stellar zero point:
    weights = rvc_err.^-2;
    if strcmp(RVtype,'serval') && ~strcmp(arm,'all')
        rvc_mean = nanwmean(rvc,weights);
        rvc = rvc - rvc_mean;
    end
    
    % calculate the std of corrected RVs:
    Nrv_corr(t) = sum(isfinite(rvc));
    bias = 1-1/4/Nrv_corr(t);
    rv_std = wstd(rvc,weights)/bias;
    rv_std_vec_corr(t) = rv_std;
    rv_med_err_corr(t) = nanmedian(rvc_err);
    
   if savedata
    
    % Write an entry to .orb file (corrected RVs, NO outliers removed):
    fprintf(fid,'STAR: %s\n',starname);
    for i=1:length(bjd)
        if ~isnan(rvc(i))
            fprintf(fid,'A %f %f %f\n',bjd(i)+mod_shift,0.001*rvc(i),0.001*rvc_err(i));
        end
    end
    fprintf(fid,'END\n');
    
    % Write an entry to all_rvs file (corrected RVs, NO outliers removed):
    for i=1:length(bjd)
        if ~isnan(rvc(i))
            fprintf(fid_all,'%14s\t%13.5f\t%9.3f\t%7.3f\t%5.3f\t%5.3f\n',starname,bjd(i)+mod_shift,rvc(i),rvc_err(i),corr_t(i),corr_err(i));
        end
    end
    
    % Write a .avc.dat file (corrected RVs, NO outliers removed):
    fid_avc = fopen([dirname 'corrected_rvs/' starname '_selfbiased.avc.dat'],'w');
    for i=1:length(bjd)
        if ~isnan(rvc(i))
            fprintf(fid_avc,'%f %f %f\n',bjd(i)+mod_shift,rvc(i),rvc_err(i));
        else
            Nnoaverage(t) = Nnoaverage(t)+1;
        end
    end
    fclose(fid_avc);
   end

end

if savedata
    fclose(fid);
    fclose(fid_all);
end

%% Find again RV-loud stars and make some display:
rv_std_median_after = nanmedian(rv_std_vec_corr(novar_rv));
if showflag
    disp('After correction:');
    fprintf(1,'Median std(RV) of RV-quiet stars: %4.2f\n',rv_std_median_after);
    fprintf(1,'std(RV) threshold for an RV_loud star: %4.1f\n',final_rv_var);
    fprintf(1,'The actually used RV-std threshold for an RV_loud star: %3.1f\n',final_rv_var);
end    

%% Plot the histogram of std(RV) after correction
if showflag
     figure(16);
     hold on
     g = logspace(0,1.05*log10(max(rv_std_vec_corr)),ceil(10*(1.05*log10(max(rv_std_vec_corr)))));
     ha = hist(rv_std_vec_corr(Nrv>=Nrv_min),g);
     stairs(g,ha,'b','LineWidth',2);
     xlim([1 10^(1.05*log10(max(rv_std_vec_corr)))]);
     ylim([0 1.05*max([ha(:);hb(:)])]);
     set(gca,'Xscale','log');
     set(gca,'xTickLabel',[1 10 100 1000 10000]);
      ylabel('Number of stars');
     legend(['Median(std(RV)) of RV-quiet stars before correction:' num2str(rv_std_median,3) ' m/s'],...
            ['Median(std(RV)) of RV-quiet stars after correction:' num2str(rv_std_median_after,3) ' m/s']);
    if saveflag
        filename = ['harps_' arm '_' RVtype '_stdRV_hist_selfbiased_' day_str];
        saveas(gcf,[dirname 'figures\' filename save_fmt]);
    end   
end

%% plot std(RV) before and after correction for RV-quiet stars:
if showflag
    figure(17)
    plot(rv_std_vec(Nrv>10 & tspan>7 & rv_std_vec<final_rv_var),rv_std_vec_corr(Nrv>10 & tspan>7 & rv_std_vec<final_rv_var),'.b','MarkerSize',15);
    hold on
    ylim([0 final_rv_var])
    xlim([0 final_rv_var])
    plot(0:1:final_rv_var,0:1:final_rv_var,'-k','LineWidth',2);
    xlabel('std(RV) before correction [m/s]');
    ylabel('std(RV) after correction [m/s]');
    grid on
    set(17,'units','normalized','outerposition',[0.25 0.05 0.70 0.95]);
    set(gca,'FontSize',16);
    title([date ': NZP correction of harps-' arm '-' RVtype ' RVs']);
    if saveflag
         filename = ['harps_' arm '_' RVtype '_RVstd_compare_selfbiased_' day_str];
         saveas(gcf,[dirname 'figures\' filename save_fmt]);
    end 
end

%% ===========================================================================
% STARTING POINT OF CORRECTING THE RADIAL VELOCITIES WITHOUT SELF BIASING
%===========================================================================
else % = if unselfbias

if savedata
    orbname = [dirname 'orb_files\harps_' arm '_' RVtype '_avc_' day_str '.orb'];
    fid = fopen(orbname,'w');
    allname = [dirname 'txt_files\harps_' arm '_' RVtype '_all_rvs_' day_str '.tsv'];
    fid_all = fopen(allname,'w');
    propname = [dirname 'txt_files\harps_' arm '_' RVtype '_empty_files_' day_str '.txt'];
    fid_prop = fopen(propname,'w');
end

rv_std_vec_corr = nan(length(targets),1);
rv_med_err_corr = nan(length(targets),1);

for t=1:target_num
    % derive the starname:
    starname = stars{t};
    fprintf(1,'\nRunning star number %d/%d (%s)\n',t,target_num,starname);
    
    % Find variable stars:
    var_rv = find(rv_std_vec>init_rv_var);
    
    % Remove also the current star:
    var_rv = [var_rv;t];
    
    % remove multiple entries:
    var_rv = sort(var_rv);
    var_rv(diff(var_rv)==0)=[];
    
%% Average RVs from the same night using all stars observed in that night:

% initialize some vectors:
nights_mean = nan(1,m);
nights_median = nights_mean;
nights_mean_err = nights_mean;
nights_mean_std = nights_mean;
Nrv_night = nan(1,m);
Nout_night = Nrv_night;

% load the results of stage 1 (the sparse RV-night matrices):
load([dirname 'mat_files\jd_mat.mat'],'jd_mat','rv_mat','err_mat');

% Remove the RVs of targets with less than Nrv_min exposures:
jd_mat(Nrv<Nrv_min,:) = nan;
rv_mat(Nrv<Nrv_min,:) = nan;
err_mat(Nrv<Nrv_min,:) = nan;

% Remove the RVs of variable targets (including the current star):
jd_mat(var_rv,:) = nan;
rv_mat(var_rv,:) = nan;
err_mat(var_rv,:) = nan;  

%% =========================================================================
% STARTING POINT OF CALCULATING THE CORRECTION (NESTED LOOP)
%===========================================================================
for y=1:m
    night_jds = jd_mat(:,y);
    night_rvs = rv_mat(:,y);
    night_err = err_mat(:,y);
    
    % remove outliers per night (but plot them) only if there are enough RVs:
    good_rvs_tmp = find(isfinite(night_rvs));
    Nrv_night_tmp = length(good_rvs_tmp);
    if Nrv_night_tmp>Nrv_min_night
        night_weights = night_err.^-2;
        rv_mean_night = nanwmean(night_rvs,night_weights);
        bias = 1-1/4/Nrv_night_tmp;
        rv_std_night = wstd(night_rvs,night_weights)/bias;
%         rv_std_night = 1.478*mad(night_rvs,1)/bias;
        bad_rv_night = find(abs(night_rvs-rv_mean_night)>sigma_outlier_night*rv_std_night & abs(night_rvs-rv_mean_night)>init_rv_var);
        night_jds(bad_rv_night) = nan;
        night_rvs(bad_rv_night) = nan;
        night_err(bad_rv_night) = nan;
        Nout_night(y) = length(bad_rv_night);
    end
    
    %==========================================================
    % Calculate the nightly zero-point RV using only good RVs:
    good_rvs = find(isfinite(night_rvs));
    Nrv_night(y) = length(good_rvs);
    bias = 1-1/4/Nrv_night(y);
    weights = night_err.^-2;
    if sum(isfinite(weights))
        nights_mean(y) = nanwmean(night_rvs(good_rvs),weights(good_rvs));
        nights_mean_err(y) = 1/sqrt(sum(weights(good_rvs)))/bias;
        nights_median(y) = wmedian(night_rvs(good_rvs),weights(good_rvs));
        nights_mean_std(y) = wstd(night_rvs(good_rvs),weights(good_rvs))/sqrt(Nrv_night(y))/bias;
        nights_mean_err(y) = max(nights_mean_err(y),nights_mean_std(y));
    end
    %==========================================================
    
end
%% ===========================================================================
% END POINT OF CALCULATING THE CORRECTION (NESTED LOOP)
%===========================================================================

%% display some statistics:
    weights_nzp = nights_mean_err.^-2;
    wmean_nights_mean = nanwmean(nights_mean(Nrv_night>=Nrv_min_night),weights_nzp(Nrv_night>=Nrv_min_night));
    std_nights_mean = wstd(nights_mean(Nrv_night>=Nrv_min_night),weights_nzp(Nrv_night>=Nrv_min_night));
    median_nights_mean_err = nanmedian(nights_mean_err(Nrv_night>=Nrv_min_night));
    fprintf(1,'std(NZP) = %4.2f m/s; med(eNZP) = %4.2f m/s\n',std_nights_mean,median_nights_mean_err);
    
%% find nights where NZP correction could not be applied and mark them:
bad_mean = find(Nrv_night<Nrv_min_night);
nights_vec = 1:m;

% Put NaN for no-mean nights for the smoothing
nights_mean(bad_mean) = NaN;
nights_mean_err(bad_mean) = NaN;

%% ==========================================================
% smooth the NZPs:
    nzp = nights_mean(:);
    enzp = nights_mean_err(:);
    weight = 1./enzp.^2;
    N_sigma_0 = abs(nzp./enzp);
    prob_n_0 = 2*(1-cdf('t',N_sigma_0',Nrv_night-1,0,1));
    [nzp_filt,nzp_filt_err] = wmeanfilt02(nights_vec,nzp,weight,smooth_span_1);
    if ~jump_corr && smooth_nzp
        if correct_out_nzp
            N_sigma = abs((nzp-nzp_filt)./enzp);
            prob_n = 2*(1-cdf('t',N_sigma',Nrv_night-1,0,1));
            outlier_nzp = find(prob_n<pvalue_outlier_nzp & prob_n_0<pvalue_outlier_nzp);
            Nout = length(outlier_nzp);
            nzp_out = nzp;
            enzp_out = enzp;
            nzp(outlier_nzp) = nan;
            enzp(outlier_nzp) = nan;
            weight = 1./enzp.^2;
            [nzp_filt,nzp_filt_err] = wmeanfilt02(nights_vec,nzp,weight,smooth_span_2);
        end
        % correct all RVs with the smoothing filter
        nights_mean = nzp_filt';
        nights_mean_err = nzp_filt_err';
        % but keep significant NZPs:
        if correct_out_nzp
            fprintf(1,'correcting for %d significant NZPs\n',Nout);
            nights_mean(outlier_nzp) = nzp_out(outlier_nzp);
            nights_mean_err(outlier_nzp) = enzp_out(outlier_nzp);
            nzp(outlier_nzp) = nzp_out(outlier_nzp);
            enzp(outlier_nzp) = enzp_out(outlier_nzp);
        end
        
    elseif jump_corr && smooth_nzp
        [nzp_filt_b,nzp_filt_err_b] = wmeanfilt02(nights_vec(1:(jd_jump-jd_min)),nzp(1:(jd_jump-jd_min)),weight(1:(jd_jump-jd_min)),smooth_span_1);
        [nzp_filt_a,nzp_filt_err_a] = wmeanfilt02(nights_vec(jd_jump-jd_min+1:end),nzp(jd_jump-jd_min+1:end),weight(jd_jump-jd_min+1:end),smooth_span_1);
        nzp_filt = [nzp_filt_b;nzp_filt_a];
        nzp_filt_err = [nzp_filt_err_b;nzp_filt_err_a];
        if correct_out_nzp
            N_sigma = abs((nzp-nzp_filt)./enzp);
            prob_n = 2*(1-cdf('t',N_sigma',Nrv_night-1,0,1));
            outlier_nzp = find(prob_n<pvalue_outlier_nzp);
            Nout = length(outlier_nzp);
            nzp_out = nzp;
            enzp_out = enzp;
            nzp(outlier_nzp) = nan;
            enzp(outlier_nzp) = nan;
            weight = 1./enzp.^2;
            [nzp_filt_b,nzp_filt_err_b] = wmeanfilt02(nights_vec(1:(jd_jump-jd_min)),nzp(1:(jd_jump-jd_min)),weight(1:(jd_jump-jd_min)),smooth_span_2);
            [nzp_filt_a,nzp_filt_err_a] = wmeanfilt02(nights_vec(jd_jump-jd_min+1:end),nzp(jd_jump-jd_min+1:end),weight(jd_jump-jd_min+1:end),smooth_span_2);
            nzp_filt = [nzp_filt_b;nzp_filt_a];
            nzp_filt_err = [nzp_filt_err_b;nzp_filt_err_a];
        end
        % correct all RVs with the smoothing filter
        nights_mean = nzp_filt';
        nights_mean_err = nzp_filt_err';
        % but keep significant NZPs:
        if correct_out_nzp
            fprintf(1,'\ncorrecting for %d significant NZPs\n',Nout);
            nights_mean(outlier_nzp) = nzp_out(outlier_nzp);
            nights_mean_err(outlier_nzp) = enzp_out(outlier_nzp);
            nzp(outlier_nzp) = nzp_out(outlier_nzp);
            enzp(outlier_nzp) = enzp_out(outlier_nzp);
        end
        
    end
        
    % optimum NZP correction: use all NZPs with error < scatter(res):
    res_rms = wstd(nzp' - nights_mean,enzp.^-2);
    ind_optim = find(enzp<res_rms);
    ind_not_optim = find(enzp>res_rms);
    if nzp_optim
        fprintf(1,'\ncorrecting for %d/%d NZPs with errors<%4.2f m/s\n',length(ind_optim),length(ind_optim)+length(ind_not_optim),res_rms);
        nights_mean(ind_optim) = nzp(ind_optim);
        nights_mean_err(ind_optim) = enzp(ind_optim);
        nights_mean_err(ind_not_optim) = res_rms;
    end
    
    % put the NZP filter as estimate and scatter as error for no-mean nights: 
        nzp_filt_est = nzp_filt(:);
        nights_mean(bad_mean) = nzp_filt_est(bad_mean);
        fprintf(1,'\ncorrecting for %d NZPs from nights with Nrv<%d by the wmean filter +/- %4.2f\n',length(bad_mean),Nrv_min_night,res_rms);
        nights_mean_err(bad_mean) = res_rms;
        
    % replace 0 and NaN NZPs (where there is no information) with the mean and rms:
        ind_nan = find(nights_mean==0 | ~isfinite(nights_mean) | ~isfinite(nights_mean_err));
        Nnan = length(ind_nan);
        if correct_nan_nzp
            fprintf(1,'\ncorrecting for %d NaN (or 0) NZPs by %4.2f +/- %4.2f\n\n',Nnan,wmean_nights_mean,std_nights_mean);
            nights_mean(ind_nan) = wmean_nights_mean;
            nights_mean_err(ind_nan) = std_nights_mean;
        end
% ==========================================================

%% calculate the global nightly drift
rv_corr = rv_mat;
for j=1:m;
    rv_corr(:,j)=rv_corr(:,j)-nights_mean(j);
end

jd = jd_mat(:);
rva = rv_corr(:);
err = err_mat(:);
inan = isnan(jd);
jd(inan)=[];
rva(inan)=[];
err(inan)=[];
clear inan
modbjd = mod(jd,1)-0.5; % jd is bjd-mod_shift; 
% =====================================================================
if drift_corr
   [a1,b1,da1,db1] = linfit(modbjd,rva,err,0);
end

% =====================================================================
%% apply the nightly correction, add its uncertainties, and write a new *numeric*, .orb, and avc.dat files:
   
% Re-load all the RVs from the dat file for correction:
    filename = targets{t};
    starname = stars{t};
    
    try
        M = load(filename);
    catch
        warning([filename ' does not exist. Press a key to continue']);
        if pauseflag; pause; end
        continue
    end
    
    % Empty files can be skipped:
    if isempty(M)
        if display_flag
            display([filename ' is an empty file. Press a key to continue']);
        end
        if pauseflag; pause; end
        continue
    end
    
    % sort the matrix:
    M = sortrows(M,1);
    bjd_tmp = M(:,1)-mod_shift;
    
    % read the data into vectors:
    bjd = M(:,1)-mod_shift; % Subtract mod_shift day from time stamps to have all
                      % the RVs of each night with the same floor(bjd).
    if strcmp(RVtype,'serval')
        rvc = M(:,6);
        rvc_err = M(:,7);
    elseif strcmp(RVtype,'drs') % drs RVs:
        rvc = M(:,8);
        rvc_err = M(:,9);
    elseif strcmp(RVtype,'crx') % crx
        rvc = M(:,12);
        rvc_err = M(:,13);
    end
    
    drift = M(:,31);
    drift(isnan(drift)) = 0;
    sa = M(:,33);
    sa(isnan(sa)) = 0;
    if given_drift_corr
       rvc = rvc - drift; 
    end
    if sa_corr
       rvc = rvc - sa; 
    end
    
    inrange = find(bjd_tmp>jd_min & bjd_tmp<jd_max);
    outrange = find(bjd_tmp<jd_min | bjd_tmp>jd_max);
        
    % Empty files should be skipped:
    if isempty(inrange)
        if display_flag
            display([filename ' has no RVs in range.']);
        end
        if pauseflag; pause; end
        continue
    end
    
    % remove values out of range:
         bjd(outrange) = [];
         rvc(outrange) = [];
         rvc_err(outrange) = [];
    
    modbjd_t = mod(bjd,1)-0.5;
    corr_t = nan(size(bjd));
    corr_err = nan(size(bjd));
    
    %==========================================================
    % subtract the nightly mean and co-add its error:
    nights = floor(bjd);
    columns = nights-jd_min+1;
    for i=1:length(bjd)
        if drift_corr
            corr_t(i) = nights_mean(columns(i)) + (b1+modbjd_t(i)*a1);
        else
            corr_t(i) = nights_mean(columns(i));
        end
        rvc(i) = rvc(i) - corr_t(i);
        if drift_corr
            corr_err(i) = sqrt(nights_mean_err(columns(i))^2 + db1^2 + (modbjd_t(i)*da1)^2);
        else
            corr_err(i) = nights_mean_err(columns(i));
        end
        rvc_err(i) = sqrt(rvc_err(i)^2 + corr_err(i)^2);
    end
    %==========================================================
    
    % for serval RVs of pre or post separately, remove again the stellar zero point:
    weights = rvc_err.^-2;
    if strcmp(RVtype,'serval') && ~strcmp(arm,'all')
        rvc_mean = nanwmean(rvc,weights);
        rvc = rvc - rvc_mean;
    end
    
    % calculate the std of corrected RVs:
    Nrv_corr(t) = sum(isfinite(rvc));
    bias = 1-1/4/Nrv_corr(t);
    rv_std = wstd(rvc,weights)/bias;
    rv_std_vec_corr(t) = rv_std;
    rv_med_err_corr(t) = nanmedian(rvc_err);
    
   % replace the nan values in M with the corrected values:
    if strcmp(RVtype,'serval')
        M(inrange,2) = rvc;
        M(inrange,3) = rvc_err;
        M(inrange,34) = corr_t;
        M(inrange,35) = corr_err;
    elseif strcmp(RVtype,'drs') % drs RVs:
        M(inrange,4) = rvc;
        M(inrange,5) = rvc_err;
        M(inrange,36) = corr_t;
        M(inrange,37) = corr_err;
    end
   
   if savedata
       
    % Write a new *numeric*.dat file
    if ~strcmp(RVtype,'crx') && ~isempty(M)
    fid_num = fopen(filename,'w');
        for i=1:size(M,1)
        fprintf(fid_num,['%13.5f %12.3f %8.3f %12.3f %8.3f  %12.3f %8.3f %12.3f %8.3f %12.3f %8.3f '...
                        '%8.3f %8.3f %8.3f %8.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %3d '...
                        '%12.3f %12.3f %12.3f %12.3f %12.3f %13.5f %10.6f %10.6f  '...
                        '%8.3f  %8.3f %8.3f  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n'],...
                        M(i,1),M(i,2),M(i,3),M(i,4),M(i,5),M(i,6),M(i,7),M(i,8),M(i,9),M(i,10),M(i,11),...
                        M(i,12),M(i,13),M(i,14),M(i,15),M(i,16),M(i,17),M(i,18),M(i,19),M(i,20),M(i,21),M(i,22),...
                        M(i,23),M(i,24),M(i,25),M(i,26),M(i,27),M(i,28),M(i,29),M(i,30),...
                        M(i,31),M(i,32),M(i,33),M(i,34),M(i,35),M(i,36),M(i,37),M(i,38),M(i,39),M(i,40),M(i,41));
        end
        fclose(fid_num);
    elseif isempty(M)
        fprintf(fid_prop,'%s\n',filename);
    end
    
    % Write an entry to .orb file (corrected RVs, NO outliers removed):
    fprintf(fid,'STAR: %s\n',starname);
    for i=1:length(bjd)
        if ~isnan(rvc(i))
            fprintf(fid,'A %f %f %f\n',bjd(i)+mod_shift,0.001*rvc(i),0.001*rvc_err(i));
        end
    end
    fprintf(fid,'END\n');
    
    % Write an entry to all_rvs file (corrected RVs, NO outliers removed):
    for i=1:length(bjd)
            fprintf(fid_all,'%14s\t%13.5f\t%9.3f\t%7.3f\t%5.3f\t%5.3f\n',starname,bjd(i)+mod_shift,rvc(i),rvc_err(i),corr_t(i),corr_err(i));
    end
    
    % Write a .avc.dat file (corrected RVs, NO outliers removed):
    fid_avc = fopen([dirname 'corrected_rvs/' starname '.avc.dat'],'w');
    for i=1:length(bjd)
            fprintf(fid_avc,'%f %f %f\n',bjd(i)+mod_shift,rvc(i),rvc_err(i));
        if isnan(rvc(i))
            Nnoaverage(t) = Nnoaverage(t)+1;
        end
    end
    fclose(fid_avc);
   end

end % end of loop on all targets

if savedata
    fclose(fid);
    fclose(fid_all);
    fclose(fid_prop);
end

if showflag
     figure(16);
     novar_rv = rv_std_vec<init_rv_var;
     rv_std_median = nanmedian(rv_std_vec(novar_rv));
     rv_std_median_after = nanmedian(rv_std_vec_corr(novar_rv));
     set(16,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
     hold off
     g = logspace(0,1.05*log10(max(rv_std_vec)),ceil(10*(1.05*log10(max(rv_std_vec)))));
     hb = hist(rv_std_vec(Nrv>=Nrv_min),g);
     stairs(g,hb,'r','LineWidth',2);
     hold on;
     xlim([1 10^(1.05*log10(max(rv_std_vec)))]);
     ylim([0 1.05*max(hb)]);
     set(gca,'Xscale','log');
     ylabel('Number of stars');
     xlabel('std(RV) [m/s]');
     title([date ': Histograms of std(RV) before and after correcting for the NZPs']);
     set(gca,'FontSize',16);
     grid on 
     g = logspace(0,1.05*log10(max(rv_std_vec_corr)),ceil(10*(1.05*log10(max(rv_std_vec_corr)))));
     ha = hist(rv_std_vec_corr(Nrv>=Nrv_min),g);
     stairs(g,ha,'b','LineWidth',2);
     xlim([1 10^(1.05*log10(max(rv_std_vec_corr)))]);
     ylim([0 1.05*max([ha(:);hb(:)])]);
     set(gca,'Xscale','log');
     set(gca,'xTickLabel',[1 10 100 1000 10000]);
      ylabel('Number of stars');
     legend(['Median(std(RV)) of RV-quiet stars before correction:' num2str(rv_std_median,3) ' m/s'],...
            ['Median(std(RV)) of RV-quiet stars after correction:' num2str(rv_std_median_after,3) ' m/s']);
    if saveflag
        filename = ['harps_' arm '_' RVtype '_stdRV_hist_' day_str];
        saveas(gcf,[dirname 'figures\' filename save_fmt]);
    end   
    
    figure(17)
    hold off
    plot(rv_std_vec(Nrv>10 & tspan>7 & rv_std_vec<final_rv_var),rv_std_vec_corr(Nrv>10 & tspan>7 & rv_std_vec<final_rv_var),'.b','MarkerSize',15);
    hold on
    ylim([0 final_rv_var])
    xlim([0 final_rv_var])
    plot(0:1:final_rv_var,0:1:final_rv_var,'-k','LineWidth',2);
    [a,b] = linfit(rv_std_vec(Nrv>10 & tspan>7 & rv_std_vec<final_rv_var),rv_std_vec_corr(Nrv>10 & tspan>7 & rv_std_vec<final_rv_var),rv_std_vec_corr(Nrv>10 & tspan>7 & rv_std_vec<final_rv_var)./sqrt(2.*Nrv(Nrv>10 & tspan>7 & rv_std_vec<final_rv_var)));
    plot(rv_std_vec(Nrv>10 & tspan>7 & rv_std_vec<final_rv_var),b+a*rv_std_vec(Nrv>10 & tspan>7 & rv_std_vec<final_rv_var),'-g','LineWidth',2);
    xlabel('std(RV) before correction [m/s]');
    ylabel('std(RV) after correction [m/s]');
    grid on
    set(17,'units','normalized','outerposition',[0.25 0.05 0.7 0.95]);
    set(gca,'FontSize',16);
    title([date ': NZP correction of harps-' arm '-' RVtype ' RVs']);
    if saveflag
         filename = ['harps_' arm '_' RVtype '_RVstd_compare_' day_str];
         saveas(gcf,[dirname 'figures\' filename save_fmt]);
    end 
    
end

end % end of unselfbias

%% ===========================================================================
% STARTING POINT OF SAVING THE RESULTS
%===========================================================================

% re-detect RV-loud stars for writing purposes
ind = find(rv_std_vec_corr>final_rv_var);

if savedata
% all stars
    fid_txt = fopen([dirname 'txt_files\harps_' arm '_' RVtype '_RVstd_' day_str '.txt'],'w');
    fprintf(fid_txt,'starname Nrv Nrv_corr RVstd_bef RVmederr_bef RVstd_aft RVmederr_aft\n');
    for i=1:length(Nrv)
        fprintf(fid_txt,'%s\t%d\t%d\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n',stars{i},Nrv(i),Nrv_corr(i),rv_std_vec(i),rv_med_err(i),rv_std_vec_corr(i),rv_med_err_corr(i));
    end
    fclose(fid_txt);

% make a txt table for RV-loud stars
    fid_act = fopen([dirname 'txt_files\rv-loud_stars_' day_str '.txt'],'w');
    fprintf(fid_act,'starname   Nrv    Nrv_corr     RVstd RVstd_err RV_mederr\n');
    for i=1:length(ind)
        fprintf(fid_act,'%s\t%d\t%d\t%7.1f\t%7.1f\t%7.1f\n',stars{ind(i)},Nrv(ind(i)),Nrv_corr(ind(i)),rv_std_vec_corr(ind(i)),rv_std_vec_corr(ind(i))/sqrt(2*Nrv(ind(i))),rv_med_err_corr(ind(i)));
    end
    fclose(fid_act);

% make a txt table for the Nightly averages
    if ~unselfbias
         fid_nav = fopen([dirname 'txt_files\harps_' arm '_' RVtype '_night_zero_' day_str '.txt'],'w');
         fprintf(fid_nav,'night   zero   err   Nrv flag\n');
         for i=1:length(nights_mean)
             fprintf(fid_nav,'%7d %7.3f %5.3f %3d %d \n',jd_min+nights_vec(i)-1,nights_mean(i),nights_mean_err(i),Nrv_night(i),flag_nzp(i));
         end
         fclose(fid_nav);
    end
    
% save the output data
    if unselfbias
       mat_file = [dirname 'mat_files\harps_' arm '_' RVtype '_stars_' day_str]; 
    else
       mat_file = [dirname 'mat_files\harps_' arm '_' RVtype '_stars_selfbiased_' day_str]; 
    end
   save(mat_file,'stars','Nspec','tspan','Nrv','Nrv_corr','rvc_mean_vec','rvc_mean_vec_err','Nnoaverage','Nout_star','Nout_night','ind_same_vec','rv_std_vec','rv_med_err','nights_mean','nights_mean_err','nights_mean_std','Nrv_night','rv_std_vec_corr','rv_med_err_corr','jd_mat','rv_mat','err_mat','nzp','enzp','flag_nzp');
end

%%
toc