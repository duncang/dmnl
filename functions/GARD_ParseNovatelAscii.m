function GARD_ParseNovatelAscii(filename)
% function GARD_ParseNovatelAscii(filename)
% Sorts the ascii log strings from 'filename' into individual files which
% can be read by dlmread.
%
% Written by Duncan Greer
% 1 March 2007
% Last Update: $Id: GARD_ParseNovatelAscii.m 1850 2008-07-14 04:52:47Z greerd $
%


disp(sprintf('Opening %s',filename));
fid = fopen(filename);

outfilename = 'bestpos.txt';
disp(sprintf('Opening %s',outfilename));
fout_bestpos = fopen(outfilename,'a');

linecount = 0;
bestposcount = 0;
while ~feof(fid)
    linedata = fgets(fid);
    linecount = linecount + 1;
    if strncmp(linedata,'#BESTPOSA',9)
       fwrite(fout_bestpos,linedata);
       bestposcount = bestposcount+1;
       disp(sprintf('BESTPOS Count = %d',bestposcount));
    end
    
end

fclose(fid);
fclose(fout_bestpos);
