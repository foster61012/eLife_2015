function [m,top_param,bottom_param] = getfits_v2(mask,junk)


top_mask = mask-circshift(mask,1);
bottom_mask = mask-circshift(mask,-1);
y = 1:length(junk);
top=[];y_top=[]; bottom=[]; y_bottom=[];

for ii = 1:length(top_mask);
    tjunk=top_mask(:,ii);
    if max(max((tjunk==max(tjunk))))==0
    else
    for jj = 1:length(y)
        if tjunk(jj) == max(tjunk)
            top=[top jj];
            y_top = [y_top y(ii)];
            break
        end
    end
    end
end
    
for ii = 1:length(bottom_mask);
    tjunk=bottom_mask(:,ii);
    if max(max(tjunk==max(tjunk)))==0
    else
    for jj = 1:length(y)
        jj = length(tjunk)-(jj-1);
        if tjunk(jj) == max(tjunk)
            bottom=[bottom jj];
            y_bottom = [y_bottom y(ii)];
            break
        end
    end
    end
end
top_param = polyfit (y_top,top,1);
bottom_param = polyfit (y_bottom,bottom,1);
% bottom_param = polyfit (y_bottom,bottom,1);
% ttt=top_param(1)*y_bottom;
%  bottom_param = mean(bottom-ttt);
 
 
m_top = top_param(1);
m_bottom = bottom_param(1);
m = mean([m_top m_bottom]);