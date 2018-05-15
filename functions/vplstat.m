function vplstat(trsid,path,VAL1,VAL2)
%*     Copyright c 1998 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Todd Walter at:      *
%*     walter@relgyro.stanford.edu                                       *
%VPLSTAT
%   VPLSTAT(FNAME, PATH, VAL1, VAL2) reads the given file of TMS Vertical 
%   Protection Limit statistics and plots the color coded histogram.  
%   In addition the regions of integrity failures and unavailability
%   are shown.  Multiple files can be combined by putting in vector of
%   strings for the FNAME (e. g. vplstat(['0421'; '0121'; '1621'])).
%   The optionalPATH is also a string and may also be in vector form 
%   (e. g. vplstat('0421',['.///////'; '../0624/']) strings must be
%   the same length). Default is './'.
%   The optional VAL1 and VAL2 arguments are given in meters (defaults
%   are 12m and 20m respectively).
%   See also : TMSSTAT, TRSNAMES, VERTSTAT, HPLSTAT
%   Example code for generating the file is included at the bottom of vplstat.m
%
%   NOTE: To plot on a non-color printer type:
%   >> colormap('gray')

%   Todd Walter May 4th 1998
%    7 July 1998 AJHansen, added a VAL argument
%   14 July 1998 Todd Walter added variable ranges
%   23 Sept 1998 AJHansen, added the total seconds variable to the input 
%    8 Oct. 1998 Todd Walter raised type above data points & added shadows
%   21 Oct. 1998 Todd Walter added vector to distinguish diagonal elements
%    3 Nov. 1998 Todd Walter added multiple input file capability
%    7 July 1999 Todd Walter changed to two VPLs, added vector paths, and
%           lines of constant probability (actual versus expected)
%   20 Oct. 1999 Todd Walter added 68%, 95%, and 99.9% stats for accuracy
%           and availability

if nargin < 1
    error('Must input a TRS id as a string or vector of strings')
end
if nargin < 2
    path = './';
end
if nargin < 3
    VAL1 = 12;
end
if nargin < 4
    VAL2 = 20;
end
[n_id id_len]=size(trsid);
[n_path path_len]=size(path);
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
  for path_idx = 1:n_path
    fname = upper(trsid(trs_idx,:));
    eval (['fid=fopen(''',path(path_idx,:),'vpl_',fname,'.hst'',''r'');']);
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
%      // AJH term to try and find/test diagonal handling
      diag_cnt = sum(diagonal);
    end
    fclose(fid);
  end
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
[i,j]    = find(data);
face_mat = [[1 2 6 5]' [2 3 7 6]' [3 4 8 7]' ...
            [4 1 5 8]' [1 2 3 4]' [5 6 7 8]']';
colors   = colormap;
for idx = 1:length(i)
  z       = log10(data(i(idx),j(idx)));
  vtx_mat = [err_bin(j(idx)) + [-0.5 0.5 0.5 -0.5 -0.5 0.5 0.5 -0.5]'*d_err_bin ...
             sig_bin(i(idx)) + [-0.5 -0.5 0.5 0.5 -0.5 -0.5 0.5 0.5]'*d_sig_bin ...
             [0 0 0 0 z z z z]'];
  c_idx   = ceil(63*(log10(data(i(idx),j(idx)))/log10(z_up_bnd))) + 1;
  patch('Vertices',  vtx_mat, ...
        'Faces',     face_mat, ...
        'FaceColor', colors(c_idx,:), ...
        'EdgeColor', 'none');
end

% determine availability and # of integrity failures
i_diag1 = find(err_bin == VAL1 | err_bin == VAL2);
i_diag2 = find(err_bin < VAL1);
i_diag3 = find(err_bin > VAL2);
i_diag4 = find(err_bin > VAL1 & err_bin < VAL2);

i_fail1 = find((err_bin(j) >= VAL1 & sig_bin(i) < VAL1) |...
               (err_bin(j) >= VAL2 & sig_bin(i) < VAL2));
n_fail1 = sum(sum(diag(data(i(i_fail1),j(i_fail1)))))...
         - sum(diagonal(i_diag1));
i_fail2 = find(err_bin(j)./sig_bin(i) >=1.0 & err_bin(j) < VAL1);
n_fail2 = sum(sum(diag(data(i(i_fail2),j(i_fail2)))))...
         - sum(diagonal(i_diag2));
i_fail3 = find(err_bin(j)./sig_bin(i) >=1.0 & sig_bin(i) > VAL2);
n_fail3 = sum(sum(diag(data(i(i_fail3),j(i_fail3)))))...
         - sum(diagonal(i_diag3));
i_fail4 = find(err_bin(j)./sig_bin(i) >=1.0 & sig_bin(i) > VAL1...
                & err_bin(j) < VAL2);
n_fail4 = sum(sum(diag(data(i(i_fail4),j(i_fail4)))))...
         - sum(diagonal(i_diag4));
i_cont  = find(sig_bin(i) >= VAL2);
n_cont  = sum(sum(diag(data(i(i_cont),j(i_cont)))));
%i_avail = find(err_bin(j) < VAL2 & sig_bin(i) < VAL2);

i_avail1 = find(err_bin(j)./sig_bin(i) < 1.0 & sig_bin(i) < VAL1);
n_avail1 = sum(sum(diag(data(i(i_avail1),j(i_avail1)))))...
         + sum(diagonal(i_diag2));

i_avail2 = find(err_bin(j)./sig_bin(i) < 1.0 & sig_bin(i) < VAL2);
n_avail2 = sum(sum(diag(data(i(i_avail2),j(i_avail2)))))...
         + sum(diagonal([i_diag2' i_diag4']));

% set the axes limits and color values
set(gca,'XLim',[x_lo_bnd x_up_bnd]);
set(gca,'YLim',[y_lo_bnd y_up_bnd]);
set(gca,'CLim',[z_lo_bnd z_up_bnd]);

% show the region of IPV operation
HT = text(0.57*(VAL2 - x_lo_bnd) + x_lo_bnd, ...
     0.95*(VAL2 - y_lo_bnd) + y_lo_bnd, ...
     0.25*sqrt(z_up_bnd*z_lo_bnd), ...
     '{\fontsize{14pt}IPV Operation}');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 1 1]);
set(HT, 'FontSize', 14);

text(0.57*(VAL2 - x_lo_bnd) + x_lo_bnd, ...
     0.95*(VAL2 - y_lo_bnd) + y_lo_bnd, ...
     0.25*sqrt(z_up_bnd*z_lo_bnd), ...
     '{\fontsize{14pt}IPV Operation}');
if n_avail2/epochs >= 0.999995
  HT = text(0.53*(VAL2 - x_lo_bnd) + x_lo_bnd, ...
            0.89*(VAL2 - y_lo_bnd) + y_lo_bnd, ...
            0.25*sqrt(z_up_bnd*z_lo_bnd), ...
            '> 99.999%');
  pos = get(HT, 'Position');
  pos = pos + [.05 -0.05 0];
  set(HT, 'Position', pos);
  set(HT, 'Color', [1 1 1]);
  set(HT, 'FontSize', 14);
  HT = text(0.53*(VAL2 - x_lo_bnd) + x_lo_bnd, ...
            0.89*(VAL2 - y_lo_bnd) + y_lo_bnd, ...
            0.25*sqrt(z_up_bnd*z_lo_bnd), ...
            '> 99.999%');
  set(HT, 'FontSize', 14);
else
  HT = text(0.57*(VAL2 - x_lo_bnd) + x_lo_bnd, ...
            0.89*(VAL2 - y_lo_bnd) + y_lo_bnd, ...
            0.25*sqrt(z_up_bnd*z_lo_bnd), ...
            [num2str(100.0*n_avail2/epochs,'%6.3f'), '%']);
  pos = get(HT, 'Position');
  pos = pos + [.05 -0.05 0];
  set(HT, 'Position', pos);
  set(HT, 'Color', [1 1 1]);
  set(HT, 'FontSize', 14);
  HT = text(0.57*(VAL2 - x_lo_bnd) + x_lo_bnd, ...
            0.89*(VAL2 - y_lo_bnd) + y_lo_bnd, ...
            0.25*sqrt(z_up_bnd*z_lo_bnd), ...
            [num2str(100.0*n_avail2/epochs,'%6.3f'), '%']);
  set(HT, 'FontSize', 14);
end

% show the region of CAT I operation
HT = text(0.45*(VAL1 - x_lo_bnd) + x_lo_bnd, ...
     0.93*(VAL1 - y_lo_bnd) + y_lo_bnd, ...
     0.25*sqrt(z_up_bnd*z_lo_bnd), ...
     '{\fontsize{14pt}CAT I Oper.}');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 1 1]);
set(HT, 'FontSize', 14);

text(0.45*(VAL1 - x_lo_bnd) + x_lo_bnd, ...
     0.93*(VAL1 - y_lo_bnd) + y_lo_bnd, ...
     0.25*sqrt(z_up_bnd*z_lo_bnd), ...
     '{\fontsize{14pt}CAT I Oper.}');
if n_avail1/epochs >= 0.999995
  HT = text(0.4*(VAL1 - x_lo_bnd) + x_lo_bnd, ...
            0.84*(VAL1 - y_lo_bnd) + y_lo_bnd, ...
            0.25*sqrt(z_up_bnd*z_lo_bnd), ...
            '> 99.999%');
  pos = get(HT, 'Position');
  pos = pos + [.05 -0.05 0];
  set(HT, 'Position', pos);
  set(HT, 'Color', [1 1 1]);
  set(HT, 'FontSize', 14);
  HT = text(0.4*(VAL1 - x_lo_bnd) + x_lo_bnd, ...
            0.84*(VAL1 - y_lo_bnd) + y_lo_bnd, ...
            0.25*sqrt(z_up_bnd*z_lo_bnd), ...
            '> 99.999%');
  set(HT, 'FontSize', 14);
else
  HT = text(0.45*(VAL1 - x_lo_bnd) + x_lo_bnd, ...
            0.84*(VAL1 - y_lo_bnd) + y_lo_bnd, ...
            0.25*sqrt(z_up_bnd*z_lo_bnd), ...
            [num2str(100.0*n_avail1/epochs,'%6.3f'), '%']);
  pos = get(HT, 'Position');
  pos = pos + [.05 -0.05 0];
  set(HT, 'Position', pos);
  set(HT, 'Color', [1 1 1]);
  set(HT, 'FontSize', 14);
  HT = text(0.45*(VAL1 - x_lo_bnd) + x_lo_bnd, ...
            0.84*(VAL1 - y_lo_bnd) + y_lo_bnd, ...
            0.25*sqrt(z_up_bnd*z_lo_bnd), ...
            [num2str(100.0*n_avail1/epochs,'%6.3f'), '%']);
  set(HT, 'FontSize', 14);
end

% outline the region of integrity failures
patch([VAL1 VAL1 VAL2 VAL2 x_up_bnd x_up_bnd], ...
      [y_lo_bnd VAL1 VAL1 VAL2 VAL2 y_lo_bnd], ...
      -[0.5 0.5 0.5 0.5 0.5 0.5], ...
      [1 0.1 0.1]);
HT = text(VAL2, ...
          0.5*(VAL2 - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['HMI']);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 0.1 0.1]);
HT = text(VAL2, ...
          0.5*(VAL2 - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['HMI']);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
HT = text(VAL2, ...
          0.4*(VAL2 - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['epochs: ', int2str(n_fail1)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 0.1 0.1]);
HT = text(VAL2, ...
          0.4*(VAL2 - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['epochs: ', int2str(n_fail1)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');

% outline the lowest region of VPL failures
patch([x_lo_bnd VAL1 VAL1], ...
      [y_lo_bnd VAL1 y_lo_bnd], ...
      -[0.5 0.5 0.5], ...
      [1 0.55 0.55]);
HT = text(0.67*(VAL1 - x_lo_bnd) + x_lo_bnd, ...
          0.35*(VAL1 - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['MI']);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 0.55 0.55]);
HT = text(0.67*(VAL1 - x_lo_bnd) + x_lo_bnd, ...
          0.35*(VAL1 - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['MI']);
set(HT,'HorizontalAlignment','Center');
set(HT, 'FontSize', 14);
HT = text(0.67*(VAL1 - x_lo_bnd) + x_lo_bnd, ...
          0.25*(VAL1 - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['epochs: ', int2str(n_fail2)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 0.55 0.55]);
HT = text(0.67*(VAL1 - x_lo_bnd) + x_lo_bnd, ...
          0.25*(VAL1 - y_lo_bnd) + y_lo_bnd, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['epochs: ', int2str(n_fail2)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');

% outline the middle region of VPL failures
patch([VAL1 VAL2 VAL2], ...
      [VAL1 VAL2 VAL1], ...
      -[0.5 0.5 0.5], ...
      [1 0.55 0.55]);
HT = text(0.67*(VAL2 - VAL1) + VAL1, ...
          0.39*(VAL2 - VAL1) + VAL1, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['MI']);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 0.55 0.55]);
HT = text(0.67*(VAL2 - VAL1) + VAL1, ...
          0.39*(VAL2 - VAL1) + VAL1, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['MI']);
set(HT,'HorizontalAlignment','Center');
set(HT, 'FontSize', 14);
HT = text(0.67*(VAL2 - VAL1) + VAL1, ...
          0.25*(VAL2 - VAL1) + VAL1, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['epochs: ', int2str(n_fail4)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 0.55 0.55]);
HT = text(0.67*(VAL2 - VAL1) + VAL1, ...
          0.25*(VAL2 - VAL1) + VAL1, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['epochs: ', int2str(n_fail4)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');

% outline the region of unavailability
patch([x_lo_bnd x_up_bnd x_up_bnd x_lo_bnd], ...
      [VAL2 VAL2 y_up_bnd y_up_bnd], ...
      -[0.5 0.5 0.5 0.5], ...
      [1 1 0.5]);
HT = text(0.50*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
          0.65*(y_up_bnd - VAL2) + VAL2, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          '{\fontsize{14pt}System Unavailable}');
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 1 0.5]);
HT = text(0.50*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
          0.65*(y_up_bnd - VAL2) + VAL2, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          '{\fontsize{14pt}System Unavailable}');
set(HT,'HorizontalAlignment','Center');
HT = text(0.50*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
          0.35*(y_up_bnd - VAL2) + VAL2, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['Alarm Epochs: ', int2str(n_cont)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 1 0.5]);
HT = text(0.50*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
          0.35*(y_up_bnd - VAL2) + VAL2, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['Alarm Epochs: ', int2str(n_cont)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');

% outline the region where integrity failures and unavailability overlap
patch([VAL2 x_up_bnd x_up_bnd], ...
      [VAL2 y_up_bnd VAL2], ...
      z_lo_bnd*-[0.45 0.45 0.45], ...
      [1 .55 0.3]);
HT = text(0.70*(x_up_bnd - VAL2) + VAL2, ...
          0.4*(y_up_bnd - VAL2) + VAL2, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['MI:']);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 0.55 0.3]);
HT = text(0.70*(x_up_bnd - VAL2) + VAL2, ...
          0.4*(y_up_bnd - VAL2) + VAL2, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          ['MI:']);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
HT = text(0.70*(x_up_bnd - VAL2) + VAL2, ...
          0.175*(y_up_bnd - VAL2) + VAL2, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          [int2str(n_fail3)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');
pos = get(HT, 'Position');
pos = pos + [.05 -0.05 0];
set(HT, 'Position', pos);
set(HT, 'Color', [1 0.55 0.3]);
HT = text(0.70*(x_up_bnd - VAL2) + VAL2, ...
          0.175*(y_up_bnd - VAL2) + VAL2, ...
          0.25*sqrt(z_up_bnd*z_lo_bnd), ...
          [int2str(n_fail3)]);
set(HT, 'FontSize', 14);
set(HT,'HorizontalAlignment','Center');


hold on;
grid off;
% make the grid visible over the patch regions (integrity, unavailability)
ytick=get(gca,'YTick');
xtick=get(gca,'XTick');
nytick=length(ytick);
nxtick=length(xtick);
for i=1:nytick
    plot3([x_lo_bnd x_up_bnd], ytick(i)*[1 1], [0.6 0.6], 'k:');
end
for i=1:nxtick
    plot3(xtick(i)*[1 1], [y_lo_bnd y_up_bnd], [0.6 0.6], 'k:');
end

% label the axes and add a title
xlabel('{\fontsize{18pt}Error (m)}');
set(get(gca,'XLabel'),'Position', ...
    [ 0.50*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
     -0.05*(y_up_bnd - y_lo_bnd) + y_lo_bnd z_lo_bnd]);
ylabel('{\fontsize{18pt}VPL_{WAAS} (m)}');
set(get(gca,'YLabel'),'Position', ...
    [-0.06*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
      0.50*(y_up_bnd - y_lo_bnd) + y_lo_bnd z_lo_bnd]);
if(n_id == 1)
  if sec_available == 1
    title (['Vertical Performance at ', ...
            trsnames(fname),' - 0x',fname, ...
            ' (',int2str(seconds),' seconds)']);
  else
    title (['Vertical Performance at ', ...
            trsnames(fname),' - 0x',fname, ...
            ' (',int2str(n_pts),' epochs)']);
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
     1.02*(y_up_bnd - y_lo_bnd) + y_lo_bnd z_lo_bnd]);

% put the sigma_v scale on the right hand y-axis
K_v_pa = 5.33;
for i = ceil(y_lo_bnd/K_v_pa):floor(y_up_bnd/K_v_pa)
  if abs(i*K_v_pa - (0.5*(y_up_bnd - y_lo_bnd) + y_lo_bnd)) > 0.05*(y_up_bnd - y_lo_bnd)
    plot3([0.99*(x_up_bnd - x_lo_bnd) + x_lo_bnd x_up_bnd], ...
          i*K_v_pa*[1 1], ...
          [0.65 0.65],'k');
    text(1.02*(x_up_bnd - x_lo_bnd) + x_lo_bnd, i*K_v_pa, 0.65, int2str(i));
  end
end
HT = text(1.03*(x_up_bnd - x_lo_bnd) + x_lo_bnd, ...
          0.5*(y_up_bnd - y_lo_bnd) + y_lo_bnd, ...
          0.65, ...
          '{\fontsize{14pt}\sigma}_V (m)');
set(HT,'HorizontalAlignment','Center');
set(HT,'Rotation',90);

plot3([x_lo_bnd VAL1], [VAL1 VAL1], [-0.5 -0.5],'k');
% put in lines of constant prob
cont_lines=0;
if cont_lines
 plot3([0 y_up_bnd/K_v_pa], [0 y_up_bnd], [z_up_bnd z_up_bnd],'k');
 plot3([0 2*y_up_bnd/K_v_pa], [0 y_up_bnd], [z_up_bnd z_up_bnd],'k');
 plot3([0 3.29*y_up_bnd/K_v_pa], [0 y_up_bnd], [z_up_bnd z_up_bnd],'k');
 bnd_68 =err_bin(bound2(0.68,data')) + d_err_bin/2;
 bnd_95 =err_bin(bound2(0.95,data')) + d_err_bin/2;
 bnd_999=err_bin(bound2(0.999,data')) + d_err_bin/2;
 plot3(bnd_68, sig_bin + d_sig_bin/2, z_up_bnd*ones(100,1), 'k')
 plot3(bnd_95, sig_bin + d_sig_bin/2, z_up_bnd*ones(100,1), 'k')
 plot3(bnd_999, sig_bin + d_sig_bin/2, z_up_bnd*ones(100,1), 'k')
end 
%save data data sig_bin err_bin diagonal

err_bnd_68 =err_bin(bound2(0.68,sum(data)')) + d_err_bin/2;
err_bnd_95 =err_bin(bound2(0.95,sum(data)')) + d_err_bin/2;
err_bnd_999 =err_bin(bound2(0.999,sum(data)')) + d_err_bin/2;
patch([err_bnd_68 + 0.01*(x_up_bnd - x_lo_bnd) ...
       err_bnd_68 ...
       err_bnd_68 - 0.01*(x_up_bnd - x_lo_bnd)], ...
      [y_lo_bnd 0.02*(y_up_bnd - y_lo_bnd) y_lo_bnd], ...
      [z_up_bnd z_up_bnd z_up_bnd],'k');   
HT=text(err_bnd_68, 0.03*(y_up_bnd - y_lo_bnd), z_up_bnd, '68%');
set(HT,'Rotation',90);
patch([err_bnd_95 + 0.01*(x_up_bnd - x_lo_bnd) ...
       err_bnd_95 ...
       err_bnd_95 - 0.01*(x_up_bnd - x_lo_bnd)], ...
      [y_lo_bnd 0.02*(y_up_bnd - y_lo_bnd) y_lo_bnd], ...
      [z_up_bnd z_up_bnd z_up_bnd],'k');   
HT=text(err_bnd_95, 0.03*(y_up_bnd - y_lo_bnd), z_up_bnd, '95%');
set(HT,'Rotation',90);
patch([err_bnd_999 + 0.01*(x_up_bnd - x_lo_bnd) ...
       err_bnd_999 ...
       err_bnd_999 - 0.01*(x_up_bnd - x_lo_bnd)], ...
      [y_lo_bnd 0.02*(y_up_bnd - y_lo_bnd) y_lo_bnd], ...
      [z_up_bnd z_up_bnd z_up_bnd],'k');   
HT=text(err_bnd_999, 0.03*(y_up_bnd - y_lo_bnd), z_up_bnd,'99.9%');
set(HT,'Rotation',90);

sig_bnd_68 =sig_bin(bound2(0.68,sum(data')', epochs)) + d_sig_bin/2;
sig_bnd_95 =sig_bin(bound2(0.95,sum(data')', epochs)) + d_sig_bin/2;
sig_bnd_999 =sig_bin(bound2(0.999,sum(data')', epochs)) + d_sig_bin/2;
patch([x_up_bnd 0.98*(x_up_bnd - x_lo_bnd) + x_lo_bnd x_up_bnd], ...
          [sig_bnd_68 + 0.01*(y_up_bnd - y_lo_bnd) ...
           sig_bnd_68 ...
           sig_bnd_68 - 0.01*(y_up_bnd - y_lo_bnd)], ...
          [z_up_bnd z_up_bnd z_up_bnd],'k');   
text(0.925*(x_up_bnd - x_lo_bnd) + x_lo_bnd, sig_bnd_68, z_up_bnd, '68%');
patch([x_up_bnd 0.98*(x_up_bnd - x_lo_bnd) + x_lo_bnd x_up_bnd], ...
          [sig_bnd_95 + 0.01*(y_up_bnd - y_lo_bnd) ...
           sig_bnd_95 ...
           sig_bnd_95 - 0.01*(y_up_bnd - y_lo_bnd)], ...
          [z_up_bnd z_up_bnd z_up_bnd],'k');   
text(0.925*(x_up_bnd - x_lo_bnd) + x_lo_bnd, sig_bnd_95, z_up_bnd, '95%');
patch([x_up_bnd 0.98*(x_up_bnd - x_lo_bnd) + x_lo_bnd x_up_bnd], ...
          [sig_bnd_999 + 0.01*(y_up_bnd - y_lo_bnd) ...
           sig_bnd_999 ...
           sig_bnd_999 - 0.01*(y_up_bnd - y_lo_bnd)], ...
          [z_up_bnd z_up_bnd z_up_bnd],'k');   
text(0.9*(x_up_bnd - x_lo_bnd) + x_lo_bnd, sig_bnd_999, z_up_bnd, '99.9%');

% put the color scale up on the right hand side
H = colorbar('vert');
set(get(H,'Ylabel'),'String','{\fontsize{14pt}Number of Points per Pixel}');
set(H,'YScale','log');
set(H,'YLim',[z_lo_bnd z_up_bnd]);
set(H,'CLim', [z_lo_bnd z_up_bnd]);

if sec_available == 1
    display (['WAAS diff for ', ...
              num2str(n_pts), ...
              ' out of ', ...
              num2str(seconds), ...
              ' seconds']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EXEMPLAR C CODE for generating vpl_????.hst files

%Initialization:
%#define HIST_BINS 100
%unsigned long vpl_stat[HIST_BINS][HIST_BINS];
%unsigned long diagonal[HIST_BINS];
%memset(vpl_stat, 0, HIST_BINS*HIST_BINS*sizeof(unsigned long));
%memset(diagonal, 0, HIST_BINS*sizeof(unsigned long));

%Accumulation:
%
%int j,k;
%double V_err, V_sig;
%
%  j = (int)floor(4.0*fabs(V_err));
%  if(j >= HIST_BINS)
%      j = HIST_BINS - 1;
%  else if(j < 0)
%      j = 0;
%  k = (int)floor(4.0*5.33*V_sig);
%  if(k >= HIST_BINS)
%      k = HIST_BINS - 1;
%  else if(k < 0)
%      k = 0;
%  vpl_stat[k][j]++;
% if((k == j) &&
%    (fabs(V_err) < 5.33*V_sig)){
%   diagonal[k]++;
% }
%
%Writing to file:
%char filename[256];
%int i,ID=0x????;
%double verr_bin_value[HIST_BINS];
%double vpl_bin_value[HIST_BINS];
%FILE *fp;
%
% for(i=0; i<HIST_BINS;i++){
%    verr_bin_value[i] = ((double)(i) + 0.5)/4.0;
%    vpl_bin_value[i]  = ((double)(i) + 0.5)/4.0;
% }
% sprintf(file_name,"vpl_%04X.hst",ID);
% fp = fopen(file_name,"wb");
% fwrite(verr_bin_value, sizeof(double), HIST_BINS,fp);
% fwrite(vpl_bin_value, sizeof(double), HIST_BINS,fp);
% fwrite(vpl_stat, sizeof(unsigned long), HIST_BINS*HIST_BINS, fp);
% fwrite(expected_number_of_WAAS_epochs, sizeof(unsigned long),1,fp);
% fwrite(diagonal, sizeof(unsigned long),HIST_BINS,fp);
% fclose(fp);



