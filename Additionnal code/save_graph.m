function save_graph(fig, FileType ,name, w, h)
%SAVE_PDF: input variables: fig that want to be saved, the width and height
% of the pdf. 
FileType = ['-d', FileType]; 

x0=0;
y0=0;
width=w;
height=h; 
set(fig, 'units', 'centimeters', 'position',[x0,y0,width,height])

set(fig, 'Units','centimeters');
screenposition = get(gcf, 'Position'); % get the figure size
set(fig,...
    'PaperUnits','centimeters', ...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize', [screenposition(3:4)]); % make the print paper size fits the figure size
print(fig, '-painters', FileType , name) 

end

