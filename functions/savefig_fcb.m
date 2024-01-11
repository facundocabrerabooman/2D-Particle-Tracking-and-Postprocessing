function []=savefig_fcb(name,x_width,y_width, format,saveflag)

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

set(gcf,'PaperSize',[x_width y_width]); %set the paper size to what you want  

set(gca,'fontsize',20) % RSI 

%%%%
if ~exist('saveflag')
 if format == 'pdf'
 saveas(gcf,[name '.' format]) 
 else 
     saveas(gcf,[name '.' format]) 
 end
end 
end

