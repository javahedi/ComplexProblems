clc
global kx ky polx poly
%clear all
%some data
%[x,y] = meshgrid(0:0.2:2,0:0.2:2);
%u = cos(x).*y;
%v = sin(x).*y;

[x,y]=meshgrid(ky,kx);
u=polx;
v=poly;


%quiver plots
figure('Position',[10 10 1000 600],'Color','w');
%hax_1 = 
subplot(1,3,1);

%left version (regular)
hq1 = quiver(x,y,u,v);
axis([-pi pi -pi pi])

%get the line position (first handle)
hkid = get(hq1,'children');
X = get(hkid(1),'XData');
Y = get(hkid(1),'YData');
%axis off;
%title('Quiver - regular ','FontSize',16);

%right version (with annotation)
%hax_2 = 
subplot(1,3,2)
cmap1 = hot(200); %colormap
cmap2 = jet(200);

for ii = 1:3:length(X)-1

    headWidth = 200 * sqrt((X(ii+1)-X(ii)).^2 + (Y(ii+1)-Y(ii)).^2); % set the headWidth, function of length of arrow
    if floor(atan2(Y(ii+1)-Y(ii),X(ii+1)-X(ii))*180/pi)>0
        
        angled = abs(floor(atan2(Y(ii+1)-Y(ii),X(ii+1)-X(ii))*180/pi)) + 1; %get the angle
        %ah = annotation('arrow','Color', cmap1(angled,:),...
        ah = annotation('arrow','Color', [0 0 1],...
        'headStyle','cback1','HeadLength',50,'HeadWidth',headWidth);
    else
        angled = abs(floor(atan2(Y(ii+1)-Y(ii),X(ii+1)-X(ii))*180/pi)) + 1; %get the angle
        %ah = annotation('arrow','Color', cmap2(angled,:),...
        ah = annotation('arrow','Color', [1 0 0],...
        'headStyle','cback1','HeadLength',50,'HeadWidth',headWidth);
        
    end
    set(ah,'parent',gca);
    set(ah,'position',[X(ii) Y(ii) X(ii+1)-X(ii) Y(ii+1)-Y(ii)]);
end
axis([-pi pi -pi pi])
%axis off;
%title('Quiver - annotations ','FontSize',16);

%linkaxes([hax_1 hax_2],'xy');

subplot(1,3,3);



for ii = 1:3:length(X)-1
    ii

    headWidth = 200 * sqrt((X(ii+1)-X(ii)).^2 + (Y(ii+1)-Y(ii)).^2); % set the headWidth, function of length of arrow
    angled = abs(floor(atan2(Y(ii+1)-Y(ii),X(ii+1)-X(ii))*180/pi)) + 1; %get the angle
    ah = annotation('arrow','Color',  cmap1(angled,:),...
        'headStyle','cback1','HeadLength',50,'HeadWidth',headWidth);
   set(ah,'parent',gca);
    set(ah,'position',[X(ii) Y(ii) X(ii+1)-X(ii) Y(ii+1)-Y(ii)]);
end
axis([-pi pi -pi pi])