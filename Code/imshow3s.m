function h=imshow3s(im,z,mag)
% IMSHOW3S displays 3D images as a 2D slice in Handle Graphics figure
%    This function displays 3D images as 2D slices in the axial 
%    orientation. The user can use the mouse-scroll to scroll through
%    different slices. The window/level can be adjusted using the left
%    mouse button. Additionally, the window/level can be reset using
%    the right mouse button.
%
%    IMSHOW3S(I) displays grayscale image I. Initialize display
%    to show the center slices.
%
%    IMSHOW3S(I,Z) displays grayscale image I. The slice number is
%    set Z. If Z=[] or out of range, then default slice number is used.
%
%    IMSHOW3S(I,Z,MAG) displays grayscale image I. The slice number is
%    set S. The displayed image is magnified by MAG percent.
%    
%
% (c) Joseph Y Cheng (jycheng@mrsrl.stanford.edu) 2011
    
% SVN info: 
%     Date:     $Date: 2012-10-12 00:37:33 -0700 (Fri, 12 Oct 2012) $
%     Revision: $Revision: 1186 $
%     Author:   $Author: jycheng $
%     Id:       $Id: imshow3s.m 1186 2012-10-12 07:37:33Z jycheng $

    ha.figure = gcf;
    clf; % clear current figure.
    
    if ((ndims(im)~=2 && ndims(im) ~= 3) || ~isreal(sum(im(:))))
        error('Input must be a real 2D/3D volume');
    end

    if (nargin < 2 || isempty(z) || z<1 || z>size(im,3))
        z = round(size(im,3)/2);        
    end
    if (nargin < 3)
        ha.mag = 100;
    else
        ha.mag = mag;
    end
    
    %% Initialize maps
    ims = im(:,:,z);
    map = [min(ims(:)) max(ims(:))];
    ha.scale_wl = abs(diff(map))/1000;
    ha.map = map;
    
    %% Setup figure and plot
    sfigure(ha.figure); hold off;
    ha.fh = imshow(im(:,:,z),ha.map,'InitialMagnification',ha.mag);  
    ha.a = gca;
    set(ha.figure,'pointer','crosshair');
    set(ha.figure,'Interruptible','off');
    imshow3_plot(ha,im,z);

    %% Setup output
    if (nargout==1)
        h = ha.figure;
    end
    
function imshow3_plot(ha,im,z)
% IMSHOW3_PLOT Function that actually plots the data.
    if qfig(ha)
        return;
    end
    %figure(ha.figure); 
    set(ha.fh,'CData',im(:,:,z));
    if ha.map(1) || ha.map(2)
        set(ha.a, 'CLim',ha.map);    
    end
    xlabel(ha.a, sprintf('XY-slice: %d',z)); 
    %drawnow;
    
    %% Re-init callback function with new values
    set(ha.figure,'WindowButtonDownFcn', {@int_callback,ha,im,z});
    set(ha.figure,'WindowScrollWheelFcn',{@int_scroll_callback,ha,im,z});
        
function int_scroll_callback(a,event,ha,im,z)
% INT_SCROLL_CALLBACK callback functions that handles mouse scrolls.
%   This allows the user to scroll through slices using the mouse
%   wheel.
    imz = size(im,3);
    dz = event.VerticalScrollCount;
    if (abs(dz)>imz)
        dz = rem(dz,imz); 
    end
    if (imz<dz+z)
        z = dz+z-imz;
    elseif (1>dz+z)
        z = dz+z+imz;
    else
        z = dz+z;
    end
        
    %% update plot
    imshow3_plot(ha,im,z);

function int_callback(obj,a,ha,im,z)
% INT_CALLBACK Main callback function that handles all functionality.
    %[imx,imy,imz] = size(im);
    if qfig(ha)
        return;
    end
    x0 = (get(ha.a,'CurrentPoint'));
    y1 = x0(1,2); x1 = x0(1,1);
   
    stype = get(obj,'SelectionType');
    if strcmp(stype,'normal')
        %% For window/level
        %dbdisp('normal selection');
        %if (0<x1 && x1<=imz && 0<y1 && y1<=imy)
            wl.x = x1; wl.y = y1;            
        %end 
        %if (exist('wl','var'))                    
            set(ha.figure,'WindowButtonMotionFcn', ...
                          {@int_motion_callback,ha,im,z,wl});           
        %end
    elseif strcmp(stype,'alt')
        %% Reset window/level
        %dbdisp(sprintf('Reset window/level: %d',z));
        ims = im(:,:,z);
        map = [min(ims(:)) max(ims(:))];
        %if (0<x1 && x1<=imz && 0<y1 && y1<=imy)
            ha.map = map;
        %end   
        imshow3_plot(ha,im,z);
    end

function int_motion_callback(a,b,ha,im,z,wl)
% INT_MOTION_CALLBACK Callback function used to set window/level.
    if qfig(ha)
        return;
    end
    temp = (get(ha.a,'CurrentPoint'));
    x1 = temp(1,1); y1 = temp(1,2);
    xdiff = x1-wl.x;
    ydiff = y1-wl.y;
    
    ha.map = modifymap(ha.map,xdiff,ydiff,ha.scale_wl);
    set(ha.figure,'WindowButtonUpFcn',{@int_release_callback,ha,im,z});
    
    %% update plots
    imshow3_plot(ha,im,z);
    
        
function int_release_callback(a,b,ha,im,z)
% INT_RELEASE_CALLBACK Callback function used to stop window/leveling.
    if qfig(ha)
        return;
    end
    set(ha.figure,'WindowButtonMotionFcn','');
    set(ha.figure,'WindowButtonUpFcn','');
    %dbdisp('realease callback');
    %ha.map
    %% update plots
    imshow3_plot(ha,im,z);
    
function nmap = modifymap(map,xdiff,ydiff,scale)
% MODIFYMAP Uses distance data to change the image map.
    [w,l] = map2wl(map);
    w = w+ydiff*scale; l = l+xdiff*scale;
    nmap = wl2map(w,l);

function [w,l] = map2wl(map)
% MAP2WL map values to window/level values.
    w = map(2)-map(1);
    l = (map(2)+map(1))/2;
    
function map = wl2map(w,l)
% WL2MAP window/level values to map values.
    map = [l-w/2 l+w/2];
    
function out = qfig(ha)
    out = ~ishandle(ha.figure) || ~ishandle(ha.a) ...
        || ha.figure~=gcf || ha.a~=gca;
