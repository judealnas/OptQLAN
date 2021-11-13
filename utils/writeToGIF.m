function writeToGIF(varargin)
%%WRITETOGIF Writes a figure to a GIF file, appending or overwriting
%
%   writeToGIF(FIG, FILE, OVERWRITE) captures the figure pointed to by
%   figure handle FIG and writes (appends) the frame to FILE if OVERWRITE
%   is true (false). If writing, the default loop count is Inf. If
%   appending, the delay time is 0.5 seconds.
%
%   writeToGIF(__,'LoopCount',loopcount) specifies the number of times to
%   loop the gif before stopping. Only applies when writing, i.e. when
%   OVERWRITE == true.
%
%   writeToGIF(__,'DelayTime', delay) specifies the time to show the
%   current frame in the gif before moving to the next frame. This option
%   only has an effect when appending, i.e. OVERWRITE == false.
p = inputParser;
p.addRequired('fig',@ishandle)
p.addRequired('file',@(x)isstring(x)||(ischar(x) && isvector(x)));
p.addRequired('b_overwrite',@islogical);
p.addParameter('LoopCount', inf, @(x)isnumeric(x) && isscalar(x) && ~mod(x,1));
p.addParameter('DelayTime', 0.5, @(x)isnumeric(x) && isscalar(x));

p.parse(varargin{:});
fig = p.Results.fig;
file = p.Results.file;
b_overwrite = p.Results.b_overwrite;
LoopCount = p.Results.LoopCount;
DelayTime = p.Results.DelayTime;

frame = getframe(fig);
im = frame2im(frame);

[imind, cm] = rgb2ind(im,256);

if b_overwrite
    imwrite(imind,cm,file,'gif','LoopCount', LoopCount, 'DelayTime', DelayTime);
else
    imwrite(imind,cm,file,'gif','WriteMode', 'append', 'DelayTime', DelayTime);
end
  

end
