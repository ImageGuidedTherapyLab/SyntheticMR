function [] = dbdisp(msg)
% DBDISP Prints a debug message.
%    DBDISP(MSG) Prints MSG to screen. The function calling DBDISP 
%    is prepended to the output.
%
% (c) Joseph Y Cheng (jycheng@mrsrl.stanford.edu) 2010
    
% SVN info: 
%     Date:     $Date: 2012-02-08 16:47:41 -0800 (Wed, 08 Feb 2012) $
%     Revision: $Revision: 1058 $
%     Author:   $Author: jycheng $
%     Id:       $Id: dbdisp.m 1058 2012-02-09 00:47:41Z jycheng $

    global JYCDebug;
    
    if JYCDebug < 1
        return;
    end
    
    %% Grab debug information
    [ST,I] = dbstack;
    
    %% Check if there is a higher level calling this function
    callname = '';
    if length(ST) > 1        
        callname = ST(2).name;
    end

    if (callname)
        disp([callname '> ' char(msg)]);
    else
        disp(msg);
    end
