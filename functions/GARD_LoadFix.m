function fix = GARD_LoadFix(filename)
% function fix = GARD_LoadFix(filename)
% Reads the fix.dat file provided by Robin Peel and selects FIXes within
% the Australian FIR.
%
% Written by Duncan Greer
%
%I
%600 Version - DAFIF data cycle 200508, build 1936, metadata FixXP700, Copyright � 2005, Robin A. Peel (robin@xsquawkbox.net).   This data is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.  You should have received a copy of the GNU General Public License along with this program ("AptNavGNULicence.txt"); if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.  This product was developed using DAFIF (the Defense Aeronautical Flight Information File), a product of the US National Imagery and Mapping Agency (NIMA). NIMA requires the following warranty statements:  (A) Under 10 U.S.C. 456, no civil action may be brought against the United States on the basis of the content of a navigational aid prepared or disseminated by either the former Defense Mapping Agency (DMA) or the National Imagery and Mapping Agency (NIMA).  (B) The DAFIF product is provided "as is," and no warranty, express or implied, including, but not limited to the implied warranties of merchantability and fitness for particular purpose or arising by statute or otherwise in law or from a course of dealing or usage in trade, is made by NIMA as to the accuracy and functioning of the product. (C): Neither NIMA nor its personnel will be liable for any claims, losses, or damages arising from or connected with the use of this product.  The user agrees to hold harmless the United States National Imagery and Mapping Agency.  The user's sole and exclusive remedy is to stop using the DAFIF product.
% 22.528056 -156.170961 00MKK


fid = fopen(filename,'r');

linedata = fgetl(fid);
linedata2 = fgetl(fid);

% scan through full file
index = 0;
fix = 0;
while ~feof(fid)
    
    linedata = fgetl(fid);
    %disp(linedata);
    
    
    if str2num(linedata(1:2)) == 99
        break;
    end
    
    if length(linedata) < 28
        continue;
    end
    
    FixName = double(linedata(24:28));
    FixLat = str2num(linedata(1:10));
    FixLon = str2num(linedata(12:22));
    
    
    
    if(FixLat > -10 || FixLat < -50)
        continue;
    end
    
    if(FixLon < 110 || FixLon > 160)
        continue;
    end
    
    index = index+1;
    fix(index,1:5) = FixName;
    fix(index,6) = FixLat*pi/180;
    fix(index,7) = FixLon*pi/180;
end





fclose(fid);
