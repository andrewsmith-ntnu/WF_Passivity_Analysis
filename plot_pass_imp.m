function []=plot_pass_imp(f,Z1,Z2,n_o,n_i,plot_type)

figure(f);

% Reset colorIdx if figure is new
persistent colorIdx colorFig
if isempty(colorFig) || ~isvalid(colorFig) || colorFig ~= gcf
    colorIdx = 1;
    colorFig = gcf;
end

t = findall(gcf, 'Type', 'tiledlayout');
if isempty(t)
    t=tiledlayout(2,1);
end



mag=abs(Z2.ResponseData);
im_part=imag(Z2.ResponseData);
pass=mimo_passivity(Z1);

%Find zero-crossings in passivity
signal = pass(:);
signChanges = diff(sign(signal));
zeroCrossings = find(signChanges ~= 0);


%Plot passivity of Z1
nexttile(t,1);
ax = gca;
colorOrder = ax.ColorOrder;
color = colorOrder(mod(colorIdx-1, size(colorOrder,1)) + 1, :);
colorIdx = colorIdx + 1;
fout=Z1.Frequency;
p=semilogx(fout,pass,'Color',color); % Use selected color
hold on
ylabel('Passivity');
grid on;
p.LineWidth=1.5;
line([Z1.Frequency(1) Z1.Frequency(end)],[0 0],'Color','k');
ax=gca;
f_size=16;
ax.FontSize=f_size;
for i = 1:length(zeroCrossings)
    xline(Z1.Frequency(zeroCrossings(i)), 'k--');  
end

%Plot magnitude/imaginary part of Z2
nexttile(t,2);
if exist('plot_type') && strcmp(plot_type,'imag')
    p=semilogx(fout,squeeze(im_part(n_o,n_i,:)),'Color',color); % Use selected color
    ylabel(['$Im[Z_{',Z2.InputName{n_i},'-',Z2.OutputName{n_o},'}]$'],'Interpreter','latex');
    hold on
else
    p=loglog(fout,squeeze(mag(n_o,n_i,:)),'Color',color); % Use selected color
    ylabel(['$|Z_{',Z2.InputName{n_i},'-',Z2.OutputName{n_o},'}|$'],'Interpreter','latex');
end


hold on;grid on

xlabel('Frequency (Hz)','interpreter','latex');
p.LineWidth=1.5;
ax=gca;
f_size=16;
ax.FontSize=f_size;


for i = 1:length(zeroCrossings)
    xline(Z1.Frequency(zeroCrossings(i)), 'k--');  
end
