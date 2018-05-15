function hplstat(trsid,HAL,path)
%*     Copyright c 1998 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Todd Walter at:      *
%*     walter@relgyro.stanford.edu                                       *
%HPLSTAT
%   HPLSTAT(FNAME, HAL, PATH) reads the given file of TMS Horizontal 
%   Protection Limit statistics and plots the color coded histogram.  
%   In addition the regions of integrity failures and unavailability
%   are shown.  The optional HAL argument is given in meters (default 30m).
%   See also : TMSSTAT, TRSNAMES, VERTSTAT, VPLSTAT
%   Example code for generating the file is included at the bottom of hplstat.m
%
%   NOTE: To plot on a non-color printer type:
%   >> colormap('gray')

%   Todd Walter May 4th 1998
%    7 July 1998 AJHansen, added a HAL argument
%   14 July 1998 Todd Walter added variable ranges
%   23 Sept 1998 AJHansen, added the total seconds variable to the input 
%    8 Oct. 1998 Todd Walter raised type above data points & added shadows
%   21 Oct. 1998 Todd Walter added vector to distinguish diagonal elements
%    3 Nov. 1998 Todd Walter added multiple input file capability


if nargin < 1
    error('Must input a TRS id as a string or vector of strings')
end
if nargin < 2
    HAL = 30;
end
if nargin < 3
    path = './';
end
[n_id id_len]=size(trsid);
%if ~ischar(path) | ~ischar(trsid) | id_len ~= 4
if ~ischar(path) | ~ischar(trsid) 
    error('Inputs must be string variables')
end

c = 1;
colors = 'brygcmw';
clf
% load in data
err_bin  = zeros(100, 1);
sig_bin  = zeros(100, 1);
data     = zeros(100, 100);
seconds  = 0;
sec_available = 1;
diagonal = zeros(100, 1);

for trs_idx =1:n_id
  fname = upper(trsid(trs_idx,:));
  eval (['fid=fopen(''',path,'hpl_',fname,'.hst'',''r'');']);
  tmp_err_bin = fread(fid,[100 1],'double');
  if isequal(err_bin,zeros(100,1))
    err_bin = tmp_err_bin;
  else
    if ~isequal(tmp_err_bin,err_bin)
      error('Data range must be the same')
    end
  end
  tmp_sig_bin = fread(fid,[100 1],'double');
  if isequal(sig_bin,zeros(100,1))
    sig_bin = tmp_sig_bin;
  else
    if ~isequal(tmp_sig_bin,sig_bin)
      error('Data range must be the same')
    end
  end
  data    = data + fread(fid,[100 100],'uint32');
  [tmp_seconds tmp_sec_available] = fread(fid,1,'uint32');
  if tmp_sec_available > 0
    seconds = seconds + tmp_seconds;
  end
  sec_available = sec_available & tmp_sec_available;
  [tmp_diagonal tmp_diagonal_available] = fread(fid,[100 1],'uint32');
  if tmp_diagonal_available > 0
    diagonal = diagonal + tmp_diagonal;
  end
  fclose(fid);
end
  data    = data';

% determine the number of points and axis ranges
n_pts = sum(sum(data));
if sec_available == 1
    epochs = seconds;
else
    epochs = n_pts;
end


d_err_bin = mean(diff(err_bin));
x_lo_bnd  = min(err_bin) - d_err_bin/2;
x_up_bnd  = max(err_bin) + d_err_bin/2;

d_sig_bin = mean(diff(sig_bin));
y_lo_bnd  = min(sig_bin) - d_sig_bin/2;
y_up_bnd  = max(sig_bin) + d_sig_bin/2;

z_lo_bnd  = 1;
z_up_bnd  = max(max(data));

% clear plot
clf;

% plot each data point as a pixel
[i,j]=find(data);
face_mat=[[1 2 6 5]' [2 3 7 6]' [3 4 8 7]' ...
          [4 1 5 8]' [1 2 3 4]' [5 6 7 8]']';
colors=colormap;
for idx = 1:length(i)
  z = log10(data(i(idx),j(idx)));
  vtx_mat = [err_bin(j(idx))+[-0.5 0.5 0.5 -0.5 -0.5 0.5 0.5 -0.5]'*d_err_bin ...
             sig_bin(i(idx))+[-0.5 -0.5 0.5 0.5 -0.5 -0.5 0.5 0.5]'*d_sig_bin ...
             [0 0 0 0 z z z z]'];
  c_idx = ceil(63*(log10(data(i(idx),j(idx)))/log10(z_up_bnd))) + 1;
  patch('Vertices',  vtx_mat, ...
        'Faces',     face_mat, ...
        'FaceColor', colors(c_idx,:), ...
        'EdgeColor', 'none');
end

% determine availability and # of integrity failures
i_diag1 = find(err_bin == HAL);
i_diag2 = find(err_bin < HAL);
i_diag3 = find(err_bin > HAL);

i_fail1 = find(err_bin(j) >= HAL & sig_bin(i) < HAL);
n_fail1 = sum(sum(diag(data(i(i_fail1),j(i_fail1)))))...
         - sum(diagonal(i_diag1));
i_fail2 = find(err_bin(j)./sig_bin(i) >=1.0 & err_bin(j) < HAL);
n_fail2 = sum(sum(diag(data(i(i_fail2),j(i_fail2)))))...
         - sum(diagonal(i_diag2));
i_fail3 = find(err_bin(j)./sig_bin(i) >=1.0 & sig_bin(i) >= HAL);
n_fail3 = sum(sum(diag(data(i(i_fail3),j(i_fail3)))))...
         - sum(diagonal(i_diag3));
i_cont  = find(sig_bin(i) >= HAL);
n_cont  = sum(sum(diag(data(i(i_cont),j(i_cont)))));
%i_avail = find(err_bin(j) < HAL & sig_bin(i) < HAL);
i_avail = find(err_bin(j)./sig_bin(i) < 1.0 & sig_bin(i) < HAL);
n_avail = sum(sum(diag(data(i(i_avail),j(i_avail)))))...
         + sum(diagonal(i_diag2));


% set the axes limits and color values
set(gca,'XLim',[x_lo_bnd x_up_bnd]);
set(gca,'YLim',[y_lo_bnd y_up_bnd]);
set(gca,'CLim',[z_lo_bnd z_up_bnd]);


% show the region of normal operation
HT=text(0.37*(HAL - x_lo_bnd) + x_lo_bnd, ...
     0.93*(HAL - y_lo_bnd) + y_lo_bnd, ...
     0.25*sqrt(z_up_bnd*z_lo_bnd), ...
     '{\fontsize{14pt}Normal Operation}');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 1 1]);
text(0.37*(HAL - x_lo_bnd) + x_lo_bnd, ...
     0.93*(HAL - y_lo_bnd) + y_lo_bnd, ...
     0.25*sqrt(z_up_bnd*z_lo_bnd), ...
     '{\fontsize{14pt}Normal Operation}');
if n_avail/epochs >= .999995
  HT = text(0.33*(HAL - x_lo_bnd) + x_lo_bnd, ...
            0.86*(HAL - y_lo_bnd) + y_lo_bnd, ...
            0.25*sqrt(z_up_bnd*z_lo_bnd), ...
            '> 99.999%');
  pos = get(HT, 'Position');
  pos = pos + [.05 -0.05 0];
  set(HT, 'Position', pos);
  set(HT, 'Color', [1 1 1]);
  set(HT, 'FontSize', 14);
  HT = text(0.33*(HAL - x_lo_bnd) + x_lo_bnd, ...
            0.86*(HAL - y_lo_bnd) + y_lo_bnd, ...
            0.25*sqrt(z_up_bnd*z_lo_bnd), ...
            '> 99.999%');
  set(HT, 'FontSize', 14);
else
  HT = text(0.37*(HAL - x_lo_bnd) + x_lo_bnd, ...
            0.86*(HAL - y_lo_bnd) + y_lo_bnd, ...
            0.25*sqrt(z_up_bnd*z_lo_bnd), ...
            [num2str(100.0*n_avail/epochs,'%6.3f'), '%']);
  pos = get(HT, 'Position');
  pos = pos + [.05 -0.05 0];
  set(HT, 'Position', pos);
  set(HT, 'Color', [1 1 1]);
  set(HT, 'FontSize', 14);
  HT = text(0.37*(HAL - x_lo_bnd) + x_lo_bnd, ...
            0.86*(HAL - y_lo_bnd) + y_lo_bnd, ...
            0.25*sqrt(z_up_bnd*z_lo_bnd), ...
            [num2str(100.0*n_avail/epochs,'%6.3f'), '%']);
  set(HT, 'FontSize', 14);
end


% outline the region of integrity failures
patch([HAL HAL x_up_bnd x_up_bnd], ...
      [y_lo_bnd HAL HAL y_lo_bnd], ...
      -[0.5 0.5 0.5 0.5], ...
      [1 0.1 0.1]);
HT = text(0.50*(x_up_bnd - HAL) + HAL, ...
          0.55*(HAL - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['HMI']);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 .1 .1]);
HT = text(0.50*(x_up_bnd - HAL) + HAL, ...
          0.55*(HAL - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['HMI']);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
HT = text(0.50*(x_up_bnd - HAL) + HAL, ...
          0.45*(HAL - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['epochs: ', int2str(n_fail1)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 .1 .1]);
HT = text(0.50*(x_up_bnd - HAL) + HAL, ...
          0.45*(HAL - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['epochs: ', int2str(n_fail1)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');

% outline the region of HPL failures
patch([x_lo_bnd HAL HAL], ...
      [y_lo_bnd HAL y_lo_bnd], ...
      -[0.5 0.5 0.5], ...
      [1 0.55 0.55]);
HT = text(0.67*(HAL - x_lo_bnd) + x_lo_bnd, ...
          0.35*(HAL - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['MI']);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 0.55 0.55]);
HT = text(0.67*(HAL - x_lo_bnd) + x_lo_bnd, ...
          0.35*(HAL - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['MI']);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
HT = text(0.67*(HAL - x_lo_bnd) + x_lo_bnd, ...
          0.25*(HAL - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['epochs: ', int2str(n_fail2)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 0.55 0.55]);
HT = text(0.67*(HAL - x_lo_bnd) + x_lo_bnd, ...
          0.25*(HAL - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['epochs: ', int2str(n_fail2)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');

% outline the region of unavailability
patch([x_lo_bnd x_up_bnd x_up_bnd x_lo_bnd], ...
      [HAL HAL y_up_bnd y_up_bnd], ...
      -[0.5 0.5 0.5 0.5], ...
      [1 1 0.5]);
HT = text(0.50*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
          0.70*(y_up_bnd - HAL) + HAL, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          '{\fontsize{14pt}System Unavailable}');
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 1 0.5]);
HT = text(0.50*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
          0.70*(y_up_bnd - HAL) + HAL, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          '{\fontsize{14pt}System Unavailable}');
set(HT,'HorizontalAlignment','Center');
HT = text(0.50*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
          0.55*(y_up_bnd - HAL) + HAL, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['Alarm Epochs: ', int2str(n_cont)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 1 0.5]);
HT = text(0.50*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
          0.55*(y_up_bnd - HAL) + HAL, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['Alarm Epochs: ', int2str(n_cont)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');

% outline the region where integrity failures and unavailability overlap
patch([HAL x_up_bnd x_up_bnd], ...
      [HAL y_up_bnd HAL], ...
      -z_lo_bnd*[0.45 0.45 0.45], ...
      [1 .55 0.3]);
HT = text(0.70*(x_up_bnd - HAL) + HAL, ...
          0.325*(y_up_bnd - HAL) + HAL, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['MI']);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 0.55 0.3]);
HT = text(0.70*(x_up_bnd - HAL) + HAL, ...
          0.325*(y_up_bnd - HAL) + HAL, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['MI']);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
HT = text(0.70*(x_up_bnd - HAL) + HAL, ...
          0.175*(y_up_bnd - HAL) + HAL, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['epochs: ', int2str(n_fail3)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 0.55 0.3]);
HT = text(0.70*(x_up_bnd - HAL) + HAL, ...
          0.175*(y_up_bnd - HAL) + HAL, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['epochs: ', int2str(n_fail3)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');

hold on;
grid off;
% make the grid visible over the patch regions (integrity, unavailability)
ytick  = get(gca,'YTick');
xtick  = get(gca,'XTick');
nytick = length(ytick);
nxtick = length(xtick);
for i = 1:nytick
    plot3([x_lo_bnd x_up_bnd], ytick(i)*[1 1], [0.6 0.6], 'k:');
end
for i = 1:nxtick
    plot3(xtick(i)*[1 1], [y_lo_bnd y_up_bnd], [0.6 0.6], 'k:');
end

% label the axes and add a title
xlabel('{\fontsize{18pt}Error (m)}');
set(get(gca,'XLabel'),'Position', ...
    [ 0.50*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
     -0.05*(y_up_bnd - y_lo_bnd) + y_lo_bnd, ...
     z_lo_bnd]);
ylabel('{\fontsize{18pt}HPL_{WAAS} (m)}');
set(get(gca,'YLabel'),'Position', ...
    [-0.06*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
      0.50*(y_up_bnd - y_lo_bnd) + y_lo_bnd z_lo_bnd]);
if(n_id == 1)
  if sec_available == 1
    title (['Horizontal Performance at ', trsnames(fname), ...
          ' - 0x',fname,' (',int2str(seconds),' seconds)']);
  else
    title (['Horizontal Performance at ', trsnames(fname), ...
          ' - 0x',fname,' (',int2str(n_pts),' epochs)']);
  end
else
  trsid2 = ['0'*ones(n_id,1) 'x'*ones(n_id,1) trsid ','*ones(n_id,1)...
            ' '*ones(n_id,1)];
  trsid3 = reshape(trsid2',[1 n_id*8]);
  if sec_available == 1
    title (['TRS: ' trsid3, ...
            ' (',int2str(seconds),' seconds)']);
  else
    title (['TRS: ' trsid3, ...
            ' (',int2str(n_pts),' epochs)']);
  end
end    
set(get(gca,'Title'),'FontSize',14);
set(get(gca,'Title'),'Position', ...
    [0.50*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
     1.02*(y_up_bnd - y_lo_bnd) + y_lo_bnd, ...
     z_lo_bnd]);


% put the sigma_H_major scale on the right hand y-axis
K_h_pa = 6.18;
for i = ceil(y_lo_bnd/K_h_pa):floor(y_up_bnd/K_h_pa)
  if abs(i*K_h_pa - (0.5*(y_up_bnd - y_lo_bnd) + y_lo_bnd)) > ...
     0.05*(y_up_bnd - y_lo_bnd)
    plot3([.99*(x_up_bnd - x_lo_bnd) + x_lo_bnd x_up_bnd], ...
          i*K_h_pa*[1 1], ...
          [0.65 0.65], ...
          'k');
    text(1.02*(x_up_bnd - x_lo_bnd) + x_lo_bnd, i*K_h_pa, 0.65, int2str(i));
  end
end
HT = text(1.04*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
          0.50*(y_up_bnd - y_lo_bnd) + y_lo_bnd, ...
          0.65, ...
          '{\fontsize{14pt}\sigma}_{H_{major}} (m)');
set(HT,'HorizontalAlignment','Center');
set(HT,'Rotation',90);

% put the color scale up on the right hand side
H = colorbar('vert');
set(get(H,'Ylabel'),'String','{\fontsize{14pt}Number of Points per Pixel}');
set(H,'YScale','log');
set(H,'YLim',[z_lo_bnd z_up_bnd]);
set(H,'CLim', [z_lo_bnd z_up_bnd]);

if sec_available == 1
    display (['WAAS differential for ', ...
              num2str(n_pts), ...
              ' out of ', ...
              num2str(seconds), ...
              ' seconds']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EXEMPLAR C CODE for generating hpl_????.hst files

%Initialization:
%#define HIST_BINS 100
%unsigned long hpl_stat[HIST_BINS][HIST_BINS];
%memset(hpl_stat, 0, HIST_BINS*HIST_BINS*sizeof(unsigned long));

%Accumulation:
%
%int j,k;
%double N_err, E_err;
%double N_cov, E_cov;
%double  NEx_cov, sig_H_major;
%
%  j = (int)floor(2.0*sqrt(N_err*N_err + E_err*E_err);
%  if(j >= HIST_BINS)
%      j = HIST_BINS - 1;
%  else if(j < 0)
%      j = 0;
% /* sig_H_major is the major axis of the horizontal error ellipse. 
%  * It is defined in Appendix J of the WAAS MOPS */
% sig_H_major =  sqrt(0.5*(E_cov + N_cov) + 
%                       sqrt(0.25*(E_cov - N_cov)*(E_cov - N_cov) + 
%                       NEx_cov*NEx_cov));
%  k = (int)floor(2.0*6.18*sig_H_major);
%  if(k >= HIST_BINS)
%      k = HIST_BINS - 1;
%  else if(k < 0)
%      k = 0;
%  hpl_stat[k][j]++;

%Writing to file:
%char filename[256];
%int i,ID=0x????;
%double herr_bin_value[HIST_BINS];
%double hpl_bin_value[HIST_BINS];
%FILE *fp;
%
% for(i=0; i<HIST_BINS;i++){
%    verr_bin_value[i] = ((double)(i) + 0.5)/2.0;
%    vpl_bin_value[i]  = ((double)(i) + 0.5)/2.0;
% }
% sprintf(file_name,"hpl_%04X.hst",ID);
% fp = fopen(file_name,"wb");
% fwrite(herr_bin_value, sizeof(double), HIST_BINS,fp);
% fwrite(hpl_bin_value, sizeof(double), HIST_BINS,fp);
% fwrite(hpl_stat, sizeof(unsigned long), HIST_BINS*HIST_BINS, fp);
% fclose(fp);


