%{ 
2020-07-03
examples: vis(SINOph); vis(SINOph,[-5 15]); vis(REC); vis(REC,[1.33 1.40]); 
    vis(REC1,REC2); vis(REC1,REC2,[scale]); colormap inferno;
LMB - move along the 3th dimension
RMB - rotate the axis
use "improfile" to get 1D cross-section
%}

function vis(REC,REC2,scale)
%input parsing - allow calls with any combination of arguments
if nargin == 1
    REC2 = [];
    scale = [];
elseif nargin == 2
    if length(REC2(:)) == length(REC(:))
        scale = [];
    else
        scale = REC2;
        REC2 = [];
    end
end

%if variable is complex -> show log(abs())
if ~isreal(REC)
    disp('complex array - displaying log(0.1+abs(REC))');
    REC = log(0.1+abs(REC));
    REC2 = log(0.1+abs(REC2));
    if isempty(scale) %its probably bad anyway
        scale = [-2 12];
    end
end

%adjust the scale automatically based on the first REC
if isempty(scale) 
    scale_FOV = round(size(REC)/4);
    scale_REC = REC(scale_FOV(1):end-scale_FOV(1),scale_FOV(2):end-scale_FOV(2),scale_FOV(3):end-scale_FOV(3));
    scale_REC = scale_REC(abs(scale_REC(:)-median(scale_REC(:))) > 0.005);
    if isempty(scale_REC)
        scale = [-0.003 0.003];
    else
        scale_REC = sort(scale_REC);
        scale = [scale_REC(round(0.01*end)) scale_REC(round(0.99*end))];
        scale = scale + [-0.01 0.01];
    end
    clearvars scale_REC
end

%initialize figure
f = figure('Units', 'Normalized');
z = round(size(REC,3)/2); %default Z crossection
if isempty(REC2)
    imagesc(REC(:,:,z),scale); axis image;
    set(gca,'YDir', 'reverse'); %'normal' - reconstructions are upside down
    pos = get(gca,'Position');
    cb = colorbar;
    set(gca,'Position',pos);
else
    ax = subplot(1,2,1,'align');
    imagesc(REC(:,:,z),scale); axis image;
    ax.YDir = 'reverse';
    
    ax = subplot(1,2,2,'align');
    imagesc(REC2(:,:,z),scale); axis image;
    ax.YDir = 'reverse';
    
    pos = get(gca,'Position');
    cb = colorbar;
    set(gca,'Position',pos);
end

udata.start = [];
udata.axisNames = ['Y'; 'X'; 'Z'];
udata.changeNames = 1;
udata.sensitivity = 1; %how fast the planes are scrolled
udata.ZOffset = round(size(REC,3)/2);
udata.ZOffset2 = 0;
udata.REC = REC; clear REC;
udata.REC2 = REC2; clear REC2;
set(f,'UserData',udata);
set(f,'WindowButtonDown',@mouseDown);
set(f,'WindowButtonUp',@mouseUp);
set(f,'WindowButtonMotionFcn', @mouseMove);
set(f,'WindowScrollWheelFcn',@Scroll);

function mouseDown(object, eventdata) 
    f = object;
    if strcmp(f.SelectionType, 'normal') %LMB - change the plane
        f.UserData.start = get(f, 'CurrentPoint');
        
    elseif strcmp(f.SelectionType, 'alt') %RMB - "rotate" the reconstruction
        f.UserData.changeNames = 1; %improves the performance by a LOT
        %make sure the axis directions make intuitive sense
        if f.UserData.axisNames(3) == 'Z' %z = -z
            f.UserData.REC = permute(f.UserData.REC,[3 1 2]);
            f.UserData.REC2 = permute(f.UserData.REC2,[3 1 2]);
            f.UserData.axisNames = circshift(f.UserData.axisNames,1);
        elseif f.UserData.axisNames(3) == 'X' %swap Z and Y axis
            f.UserData.REC = permute(f.UserData.REC,[1 3 2]);
            f.UserData.REC2 = permute(f.UserData.REC2,[1 3 2]);
            f.UserData.axisNames = circshift(f.UserData.axisNames,1);
            f.UserData.axisNames(1:2) = circshift(f.UserData.axisNames(1:2),1);
        elseif f.UserData.axisNames(3) == 'Y' %revert z and x case changes
            f.UserData.REC = permute(f.UserData.REC,[3 2 1]);
            f.UserData.REC2 = permute(f.UserData.REC2,[3 2 1]);
            f.UserData.axisNames(1:2) = circshift(f.UserData.axisNames(1:2),1);
            f.UserData.axisNames = circshift(f.UserData.axisNames,1);
        end
    end
    f.UserData.sensitivity = (20+size(f.UserData.REC,3))/1;
end

function mouseUp(object, eventdata)
    f = object;
    f.UserData.start = [];
    %save the offset to start scrolling from the plane when the mouse was released
    try
        f.UserData.ZOffset = f.UserData.ZOffset2;
    catch error
    end
end

function mouseMove (object, eventdata)
    f = object;
    if ~isempty(f.UserData.start)
        mousePos = get(f, 'CurrentPoint');
        shift = f.UserData.start(2) - mousePos(2); %diff [ranging from 0 to 1] in figure coordinates
        
        %start scrolling from the middle plane and clear the offset
        if f.UserData.changeNames == 1
            f.UserData.ZOffset = round(size(f.UserData.REC,3)/2);
            f.UserData.ZOffset2 = 0;
        end

        %scale to Z values
        Zindex = round(f.UserData.ZOffset + shift*f.UserData.sensitivity);
        %update the last Z index
        f.UserData.ZOffset2 = Zindex;
        %min/max boundary check
        if Zindex >= size(f.UserData.REC,3); Zindex = size(f.UserData.REC,3); end
        if Zindex < 1; Zindex = 1; end
        
        %update figure, single input case
        if isempty(f.UserData.REC2)
            ax = f.CurrentAxes;
            target = get(ax,'Children');
            set(target,'CData',f.UserData.REC(:,:,Zindex));
            set(get(ax, 'title'),'string',[f.UserData.axisNames(3) ' = ' num2str(Zindex)]);
            if f.UserData.changeNames == 1
                ax.YLim = [0.5 size(f.UserData.REC,1)];
                ax.XLim = [0.5 size(f.UserData.REC,2)];
                set(get(ax, 'ylabel'),'string',f.UserData.axisNames(1));
                set(get(ax, 'xlabel'),'string',f.UserData.axisNames(2));
                f.UserData.changeNames = 0;
            end
        else %comparison case
            axes = findall(f,'type','axes');
            axLeft = axes(2);
            target1 = get(axLeft,'Children');
            set(target1,'CData',f.UserData.REC(:,:,Zindex));
            set(get(axLeft, 'title'),'string',[f.UserData.axisNames(3) ' = ' num2str(Zindex)]);
            if f.UserData.changeNames == 1
                axLeft.YLim = [0.5 size(f.UserData.REC,1)];
                axLeft.XLim = [0.5 size(f.UserData.REC,2)];
                set(get(axLeft, 'ylabel'),'string',f.UserData.axisNames(1));
                set(get(axLeft, 'xlabel'),'string',f.UserData.axisNames(2));
            end

            axRight = axes(1);
            target2 = get(axRight,'Children');
            set(target2,'CData',f.UserData.REC2(:,:,Zindex));
            set(get(axRight, 'title'),'string',[f.UserData.axisNames(3) ' = ' num2str(Zindex)]);
            if f.UserData.changeNames == 1
                axRight.YLim = [0.5 size(f.UserData.REC,1)];
                axRight.XLim = [0.5 size(f.UserData.REC,2)];
                set(get(axRight, 'ylabel'),'string',f.UserData.axisNames(1));
                set(get(axRight, 'xlabel'),'string',f.UserData.axisNames(2));
                f.UserData.changeNames = 0;
            end
        end
    end
    pause(0.01);
end

%mouse scroll adjusts the zoom of the slice
function Scroll(object, eventdata)
    f = object;

    sensitivity = 12;
    
    y_limits = size(f.UserData.REC,1);
    x_limits = size(f.UserData.REC,2);
    mousePos = get(f, 'CurrentPoint'); %cursor location
    y_weight = (0.5-mousePos(2));%speed up the zoom in the direction which is closer to the edge of the figure
    x_weight = (0.5-mousePos(1));
    rectangularSensitivity = y_limits/x_limits;
    ScrollCount = eventdata.VerticalScrollCount * sensitivity;

    axes = findall(f,'type','axes');
    ax = axes(1);
    try
    ax.YLim = [ax.YLim(1) - (0.5+y_weight)*ScrollCount*rectangularSensitivity ax.YLim(2) + (0.5-y_weight)*ScrollCount*rectangularSensitivity];
    ax.XLim = [ax.XLim(1) - (0.5-x_weight)*ScrollCount ax.XLim(2) + (0.5+x_weight)*ScrollCount];
    catch
    end
    %check the limits
    if ax.YLim(1) < 1; ax.YLim(1) = 0.5; end
    if ax.YLim(2) > y_limits; ax.YLim(2) = y_limits+0.5; end
    if ax.XLim(1) < 1; ax.XLim(1) = 0.5; end
    if ax.XLim(2) > x_limits; ax.XLim(2) = x_limits+0.5; end

    %case with 2 images
    if ~isempty(f.UserData.REC2)
        ax = axes(2);
        try
        ax.YLim = [ax.YLim(1) - (0.5+y_weight)*ScrollCount*rectangularSensitivity ax.YLim(2) + (0.5-y_weight)*ScrollCount*rectangularSensitivity];
        ax.XLim = [ax.XLim(1) - (0.5-x_weight)*ScrollCount ax.XLim(2) + (0.5+x_weight)*ScrollCount];
        catch
        end
        if ax.YLim(1) < 1; ax.YLim(1) = 0.5; end
        if ax.YLim(2) > y_limits; ax.YLim(2) = y_limits+0.5; end
        if ax.XLim(1) < 1; ax.XLim(1) = 0.5; end
        if ax.XLim(2) > x_limits; ax.XLim(2) = x_limits+0.5; end
    end
end
        
end