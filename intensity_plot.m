% % % 
% % % Code to measure the intensity distribution of a contracting network out of a sequence of images.
% % % The output is a variable called Profile that will have a temporal sequence of the density profiles.
% % % The general approach is to fit Gaussians to determine which parts of the image are the contracting network,
% % % which parts are areas inside the channel and outside the network, and which parts are outside of the channel. 
% % % A mask is made that incorparates the network, and the intensity contribution from unpolymerized tubulin is 
% % % approximated by the peak intensity of the parts of the channel outside of the network and is used to 
% % % correct the network intensity. From there, straight lines are fit to the top and bottom of the network that are
% % % used to rotate the image, such that the intensity can be averaged along the length of the channel. The intensity
% % % profile is then normalized such that the area under the profile is equal to 1.
% % % 
% % % 


close all
directory='/Users/peterfoster/Desktop/Intensity/1.5mm/8_14_8_5extract_slide2/'; %this should be a directory containing a sequence of tif images and the auxilary matlab functions

files=dir(fullfile(directory,'*.tif'));

area_in_pixels = 225; % Image should be square, and this is length of one of the square's size (in pixels)
Means(length(files),:)=[.1 ;3 ;8];
Spreads(length(files),:)=[.03;.2; 3];
Amps(length(files),:)=[.7 ;.2;.3];
precontraction=0;


%% make mask
for i = 1:length(files)

    ii = length(files)-(i-1);
    i=ii;
i

   temp_frames(i,:,:)=importdata(strcat(directory,files(i).name));
   junk=double(squeeze(temp_frames(i,:,:)));
   junk=junk+min(min(junk));
   junk=10*junk;
   Max_int = max(max(junk));
   
   mask = zeros(size(junk));
   param = fspecial('gaussian',size(junk),1.0);


   
  if i ==length(files) 
 [count,hist_x] = hist(reshape(junk,1,225*225),1000);
 count = count/(mean(count)*length(count));
  f = fit( hist_x', count', 'a*exp(-((x-b).^2)./c) + d*exp(-((x-f).^2)./g) +  h*exp(-((x-k).^2)/l)' ,'Startpoint',[0.7 0.1 0.03 0.2 3 0.2 0.3 8 3],'Lower',[0 0 0 0 0 0 0 0 0],'Upper',[1 max(hist_x) 20 1 max(hist_x) 20 1 max(hist_x) 20],'TolFun',10^-10,'TolX',10^-10,'MaxFunEvals',500000000000);
 coef = coeffvalues(f);
means = [coef(2) coef(5) coef(8)];
amplitudes = [coef(1) coef(4) coef(7)];
spread = sqrt([coef(3) coef(6) coef(9)]);
Means(i,:)=[coef(2); coef(5); coef(8)];
Spreads(i,:) = sqrt([coef(3); coef(6); coef(9)]);
Amps(i,:) = [coef(1) ;coef(4) ;coef(7)];
[means,ind]=sort(means);
spread = spread(ind);
amplitudes=amplitudes(ind);
Background_signal=means(2)


  else
      
      
      
      
      
 [count,hist_x] = hist(reshape(junk,1,225*225),1000);
 count = count/(mean(count)*length(count));
 f = fit( hist_x', count', 'a*exp(-((x-b).^2)./c) + d*exp(-((x-f).^2)./g) +  h*exp(-((x-k).^2)/l)' ,'Startpoint',[Amps(i+1,1) Means(i+1,1) Spreads(i+1,1) Amps(i+1,2) Means(i+1,2) Spreads(i+1,2) Amps(i+1,3) Means(i+1,3) Spreads(i+1,3)],'Lower',[0 0 0 0 0 0 0 0 0],'Upper',[1 max(hist_x) 20 1 max(hist_x) 20 1 max(hist_x) 20],'TolFun',10^-10,'TolX',10^-10,'MaxFunEvals',500000000000);
 coef = coeffvalues(f);
means = [coef(2) coef(5) coef(8)];
amplitudes = [coef(1) coef(4) coef(7)];
spread = sqrt([coef(3) coef(6) coef(9)]);
Means(i,:)=[coef(2); coef(5); coef(8)];
Spreads(i,:) = sqrt([coef(3); coef(6); coef(9)]);
Amps(i,:) = [coef(1) ;coef(4) ;coef(7)];
[means,ind]=sort(means);
spread = spread(ind);
amplitudes=amplitudes(ind);
 Means_plus_std = means + spread;
 Means_minus_std = means - spread;

if  Means_plus_std(1)<Means_minus_std(2) && Means_plus_std(2) <  Means_minus_std(3) && precontraction ==0
Background_signal=means(2)   


else
    precontraction = 1;
         f = fit( hist_x', count', 'a*exp(-((x-b).^2)./c) + d*exp(-((x-f).^2)./g)','Startpoint',[Amps(i+1,1) Means(i+1,1) Spreads(i+1,1) Amps(i+1,2) Means(i+1,2) Spreads(i+1,2)],'Lower',[0 0 0 0 0 0],'Upper',[1 max(hist_x) 20 1 max(hist_x) 20],'TolFun',10^-10,'TolX',10^-10,'MaxFunEvals',5000000000000);

coef = coeffvalues(f);
means = [coef(2) coef(5)];
amplitudes = [coef(1) coef(4)];
spread = [coef(3) coef(6)];
Means(i,:)=[coef(2); coef(5);0];
Spreads(i,:) = [coef(3); coef(6);0];
Amps(i,:) = [coef(1) ;coef(4) ;0];
[means,ind]=sort(means);
spread = spread(ind);
amplitudes=amplitudes(ind);



end


end
  
      




threshold_list(i) = max(means) - 0.75*spread(means==max(means));



   mask(junk>threshold_list(i)) = 1;


   mask=imfill(mask);
   mask=bwareaopen(mask,1000);


   
   Mask(i,:,:)=mask;

%% Fit straight lines to the top and bottom of the network and find the average slope 
[m,top_param,bottom_param] = getfits_v2(mask,junk);

if (m < 0)
    junk=  fliplr(junk);
    mask = fliplr(mask);
    [m,top_param,bottom_param] = getfits_v2(mask,junk);

end




Mask1=mask;
mask2=nan(size(mask));
mask2(mask==1)=1;
mask=mask2;


%% Check average slope, and if it's too large, rotate the image. 
sum = zeros(length(1:length(junk)-1),area_in_pixels);
low_slope(i) = abs(m) < 1/area_in_pixels;
if abs(m) < 1/area_in_pixels
    d = nanmean(junk.*mask,2);
    Find_mid_list = sort(Means(i,:));
%     d = d - Find_mid_list(2);
    d(1:9) = 0;
    d(length(d)-9:end)=0;
        d=d-Background_signal;

    d=d(10:length(d)-10);
    d = d./trapz(d(~isnan(d)));  
    plot(d)
    Profile(i,:) = d;
    axis([5 area_in_pixels-5 0 .1])
    drawnow
    F(i)=getframe;
    clear d
    
elseif abs(m)>=1/area_in_pixels
for jj = 10:area_in_pixels-10
    jj=jj+j_offset;
for ii = 1:length(junk)-1;
    %initialize
  x_left_mid=0;y_left_mid=0;x_right_mid=0;y_right_mid=0;x_upper_mid=0;y_upper_mid=0;x_lower_mid=0;y_lower_mid=0; %sum=0
    
    lower_y_bound = ii;
    upper_y_bound = ii+1;
    x_left_upper_bound = m*lower_y_bound+jj;
    x_left_lower_bound = m*lower_y_bound+jj+1;
    x_right_upper_bound = m*(upper_y_bound)+jj;
    x_right_lower_bound =m*(upper_y_bound)+jj+1;
    
    %get intersection verticies
    if floor(x_left_lower_bound) - floor(x_left_upper_bound) ~=0
        x_left_mid = floor(x_left_lower_bound);
        y_left_mid = lower_y_bound; 
    end
    if floor(x_right_lower_bound) - floor(x_right_upper_bound) ~=0
        x_right_mid = floor(x_right_lower_bound);
        y_right_mid = upper_y_bound;
    end
    if floor (x_right_upper_bound) - floor(x_left_upper_bound) ~=0
        x_upper_mid = floor (x_right_upper_bound);
        y_upper_mid = (x_upper_mid - jj)/m;
        
    end
    if floor (x_right_lower_bound) - floor(x_left_lower_bound) ~=0
        x_lower_mid = floor (x_right_lower_bound);
        y_lower_mid = (x_lower_mid - jj-1)/m;
    end
    
    
    
    
    
    if abs(ceil(x_left_lower_bound) - ceil(x_left_upper_bound)) ==1
        
    %upper triangle
    if y_left_mid ~=0 && y_upper_mid~=0
          sum(ii,jj)=sum(ii,jj)+0.5*abs(y_upper_mid-lower_y_bound)*abs(x_left_mid  -  x_left_upper_bound)*junk(floor(x_left_upper_bound),lower_y_bound).*mask(floor(x_left_upper_bound),lower_y_bound);

    end
    %lower triangle
    if y_right_mid ~=0 && y_lower_mid~=0
          sum(ii,jj)=sum(ii,jj)+0.5*abs(x_right_lower_bound-x_right_mid)*junk(floor(x_right_lower_bound),lower_y_bound).*mask(floor(x_right_lower_bound),lower_y_bound);

    end
    
    %mid parallelagram if there's 2 triangles
    if y_right_mid ~=0 && y_lower_mid~=0 && y_left_mid ~=0 && y_upper_mid~=0
        area = (1-0.5*abs(upper_y_bound-y_upper_mid)*abs(x_right_upper_bound-floor(x_right_upper_bound))- 0.5*abs(y_lower_mid-lower_y_bound)*abs(ceil(x_left_lower_bound)-x_left_lower_bound));
          sum(ii,jj) = sum(ii,jj) + junk( floor(x_left_lower_bound),lower_y_bound)*area.*mask( floor(x_left_lower_bound),lower_y_bound);
    end
    
    
    
    
    if y_upper_mid == 0 && y_lower_mid==0 &&y_right_mid ~=0 && y_left_mid ~=0
        %upper rectangle
      sum(ii,jj) = sum(ii,jj) + (ceil(x_right_upper_bound) -  max(x_right_upper_bound,x_left_upper_bound)  ) *junk(floor(x_right_upper_bound),lower_y_bound).*mask(floor(x_right_upper_bound),lower_y_bound);
        %%upper triangle
        sum(ii,jj) = sum(ii,jj) + 0.5*abs(x_right_upper_bound-x_left_upper_bound)*junk(floor(x_right_upper_bound),lower_y_bound);
        % lower rectangle
        sum(ii,jj) = sum(ii,jj) +( min(x_right_lower_bound,x_left_lower_bound) - floor(x_left_lower_bound) ) *junk(floor(x_right_lower_bound),lower_y_bound).*mask(floor(x_right_lower_bound),lower_y_bound);
        % lower triangle 
        sum(ii,jj) = sum(ii,jj) + 0.5*abs(x_right_lower_bound-x_left_lower_bound)*junk(floor(x_right_lower_bound),lower_y_bound);
    end  
    
    else
        sum(ii,jj)=nan;
    end
    
   
  

        
end

end

d = nanmean(sum);
nonzero=d(d~=0);
Find_mid_list = sort(Means(i,:));
    d=d(10:length(d)-10);
        d=d-Background_signal;

d = d./trapz(d(~isnan(d))); %normalize profile such that the area underneath is set equal to 1.
plot(d)
Profile(i,:) = d;
clear d
axis([5 area_in_pixels-5 0 .05])
drawnow
F(i)=getframe;

end



end
