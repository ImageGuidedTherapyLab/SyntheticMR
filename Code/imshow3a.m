function h=imshow3(im,x,y,z)
%IMSHOW3 Displays 3D images as 2D images in Handle Graphics figure.
%   This function displays 3D images as 2D slices in the saggital,
%   coronal, and axial orientation. The user can select the slices
%   to view. Additionally, the window/level can be set using the
%   mouse.
%
%   IMSHOW3(I) displays grayscale image I. Initialize display
%   to show the center slices.
%
%   IMSHOW3(I,x,y,z) displays grayscale image I. The slices are
%   set by x, y, and z. 
%    
%
%   (c) jycheng@mrsrl.stanford.edu 2010.
    
% SVN info: 
%     Date:     $Date: 2013-12-02 16:36:02 -0800 (Mon, 02 Dec 2013) $
%     Revision: $Revision: 1380 $
%     Author:   $Author: jycheng $
%     Id:       $Id: imshow3.m 1380 2013-12-03 00:36:02Z jycheng $

    ha.figure = gcf;
    clf; % clear current figure.
    
    if (nargin < 2)
        temp = round(size(im)/2);        
        x = temp(1); y = temp(2); z = temp(3);
    end
    
    %% Initialize maps
    map = [min(im(:)) max(im(:))];
    ha.scale_wl = max(im(:))/100;
    ha.map1 = map; ha.map2 = map; ha.map3 = map;
    
    [nx,ny,nz] = size(im);
    prp.tot = ny+2*nz;
    prp.yz = [1:nz];
    prp.xz = [nz+(1:nz)];
    prp.xy = [nz*2+(1:ny)];
    ha.prp = prp;
    
    %% Create the button group.
    hb = uibuttongroup('visible','off','Position',[0 0 0.1 1],...
                       'bordertype','none',...
                       'backgroundcolor',[0.8 0.8 0.8]);
    %% Create three radio buttons in the button group.
    ha.u0 = uicontrol('Style','Radio','String','Reset win/level',...
                   'pos',[0  0 120 30],'parent',hb,'HandleVisibility','off');
    ha.u1 = uicontrol('Style','Radio','String','Window/level',...
                   'pos',[0 30 120 30],'parent',hb,'HandleVisibility','off');
    ha.u2 = uicontrol('Style','Radio','String','Slice select',...
                   'pos',[0 60 120 30],'parent',hb,'HandleVisibility','off');
    %% Initialize some button group properties. 
    set(hb,'SelectionChangeFcn',{@int_radio_callback,ha});
    set(hb,'SelectedObject',ha.u2);  
    set(hb,'Visible','on');
        
    %% Setup figure and plot
    figure(ha.figure);
    imshow3_plot(ha,im,x,y,z);
    set(ha.figure,'pointer','crosshair');

    %% Setup output
    if (nargout==1)
        h = ha.figure;
    end
    
function imshow3_plot(ha,im,x,y,z)
% IMSHOW3_PLOT Function that actually plots the data.
    
    figure(ha.figure);
    p = ha.prp;
    
    %% Saggital
    subplot(1,p.tot,p.yz); imshow(squeeze(im(x,:,:)),ha.map1); 
    title(sprintf('YZ:%d',x)); ha.a1 = gca;

    %% Coronal
    subplot(1,p.tot,p.xz); imshow(squeeze(im(:,y,:)),ha.map2); 
    title(sprintf('XZ:%d',y)); ha.a2 = gca;
    
    %% Axial
    subplot(1,p.tot,p.xy); imshow(im(:,:,z),ha.map3); 
    title(sprintf('XY:%d',z)); ha.a3 = gca;
    
    %% Re-init callback function with new values
    set(ha.figure,'WindowButtonDownFcn', {@int_callback,ha,im,x,y,z});
        
function int_callback(obj,event,ha,im,x,y,z)
% INT_CALLBACK Main callback function that handles all functionality.
    [imx,imy,imz] = size(im);
    
    x1 = round(get(ha.a1,'CurrentPoint'));
    y1 = x1(1,2); x1 = x1(1,1);
    x2 = round(get(ha.a2,'CurrentPoint'));
    y2 = x2(1,2); x2 = x2(1,1);
    x3 = round(get(ha.a3,'CurrentPoint'));
    y3 = x3(1,2); x3 = x3(1,1);
    
    if (get(ha.u1,'Value'))
        %% For window/level
        if (0<x1 && x1<=imz && 0<y1 && y1<=imy)
            wl.axis = 0; wl.x = x1; wl.y = y1;
        elseif (0<x2 && x2<=imz && 0<y2 && y2<=imx)
            wl.axis = 1; wl.x = x2; wl.y = y2;
        elseif (0<x3 && x3<=imx && 0<y3 && y3<=imy)
            wl.axis = 2; wl.x = x3; wl.y = y3;
        end 
        if (exist('wl'))
            set(ha.figure,'WindowButtonMotionFcn', ...
                          {@int_motion_callback,ha,im,x,y,z,wl});
            set(ha.figure,'WindowButtonUpFcn',{@int_release_callback,ha});
        end
    elseif (get(ha.u2,'Value'))  
        %% For slice selection        
        if (0<x1 && x1<=imz && 0<y1 && y1<=imy)
            imshow3_plot(ha,im,x,y1,x1);
        elseif (0<x2 && x2<=imz && 0<y2 && y2<=imx)
            imshow3_plot(ha,im,y2,y,x2);
        elseif (0<x3 && x3<=imx && 0<y3 && y3<=imy)
            imshow3_plot(ha,im,y3,x3,z);
        end   
    else
        map = [min(im(:)) max(im(:))/2];
        if (0<x1 && x1<=imz && 0<y1 && y1<=imy)
            ha.map1 = map;
        elseif (0<x2 && x2<=imz && 0<y2 && y2<=imx)
            ha.map2 = map;
        elseif (0<x3 && x3<=imx && 0<y3 && y3<=imy)
            ha.map3 = map;
        end   
        imshow3_plot(ha,im,x,y,z);
    end

function int_motion_callback(obj,event,ha,im,x,y,z,wl)
% INT_MOTION_CALLBACK Callback function used to set window/level.
    if (wl.axis == 0)
        temp = round(get(ha.a1,'CurrentPoint'));
        x1 = temp(1,1); y1 = temp(1,2);
        xdiff = x1-wl.x; ydiff = y1-wl.y;
        ha.map1 = modifymap(ha.map1,xdiff,ydiff,ha.scale_wl);
    elseif (wl.axis == 1)
        temp = round(get(ha.a2,'CurrentPoint'));
        x1 = temp(1,1); y1 = temp(1,2);
        xdiff = x1-wl.x; ydiff = y1-wl.y;
        ha.map2 = modifymap(ha.map2,xdiff,ydiff,ha.scale_wl);
    else
        temp = round(get(ha.a3,'CurrentPoint'));
        x1 = temp(1,1); y1 = temp(1,2); 
        xdiff = x1-wl.x; ydiff = y1-wl.y;
        ha.map3 = modifymap(ha.map3,xdiff,ydiff,ha.scale_wl);
    end

    %% update plots
    imshow3_plot(ha,im,x,y,z);
        
function int_release_callback(obj,event,ha)
% INT_RELEASE_CALLBACK Callback function used to stop window/leveling.
    set(ha.figure,'WindowButtonMotionFcn','');
    
function int_radio_callback(obj,event,ha)
% INT_RADIO_CALLBACK Callback function used to switch cursor type.    
    selectobj = get(obj,'SelectedObject');
    if (ha.u0 == selectobj)
        set(ha.figure,'pointer','arrow');
    elseif(ha.u1 == selectobj)
        set(ha.figure,'pointer','arrow');
    else
        set(ha.figure,'pointer','crosshair');
    end
    
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