function [out,out_err] = wmeanfilt02(tin,in,weight,win)
% [out,out_err] = wmeanfilt02(tin,in,weight,win)
% Moving wmean filter that uses time units,
% and is more accurate at the edges.
% 
% INPUT:
% tin - timestamps of input data
% in - input data (can contain NaN)
% weight - input weights (can contain NaN)
% win - window for calculating wmean (in tin units).
% 
% OUPUT:
% out - output wmean filter
% out_err - output wmean filter uncertainty
% 
% Modification history:
% 20180218 LT: created
% 20180430 LT: added the uncertainty as output

hwin = floor(win/2);
len = length(in);
out = nan(size(in));
out_err = out;
% weights of NaN elements must be NaN:
weight(isnan(in)) = nan;
% filter (take care of edges):
 for i=1:len;  
     if tin(i)-tin(1)<hwin
        out(i) = nanwmean(in(tin<tin(i)+hwin),weight(tin<tin(i)+hwin));
        out_err(i) = 1/sqrt(nansum(weight(tin<tin(i)+hwin)));
     elseif tin(end)-tin(i)<hwin
        out(i) = nanwmean(in(tin>tin(i)-hwin),weight(tin>tin(i)-hwin)); 
        out_err(i) = 1/sqrt(nansum(weight(tin>tin(i)-hwin)));
     else
        out(i) = nanwmean(in(tin>tin(i)-hwin & tin<tin(i)+hwin),weight(tin>tin(i)-hwin & tin<tin(i)+hwin));
        out_err(i) = 1/sqrt(nansum(weight(tin>tin(i)-hwin & tin<tin(i)+hwin)));
     end
 end
