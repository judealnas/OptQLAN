function writeToGIF(varargin)
%%WRITETOGIF writes a figure to a GIF file, appending or overwriting
p = inputParser;
p.addRequired('fig',@ishandle)
p.addRequired('file',@isstring);
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
