function y=participation_factors_VSC_GF(conv,conv_str)
close all
pole_number=size(conv.A,1);

%% calculate right and left eigenvalues
[R,D]=eig(conv.A);
L=inv(R);
eig_A=diag(D);
figure(1)
plot(real(eig_A),imag(eig_A),'ro')
hold on
grid on

% calculate participation factors
participation_factor=L.*R.';

for i_pole=1:pole_number
    
    % if real(eig_A(i_pole))>0
    figure(2)
    subplot(2,1,1)
    plot(real(eig_A),imag(eig_A),'o')
    hold on
    if real(eig_A(i_pole))>0
        plot(real(eig_A(i_pole)),imag(eig_A(i_pole)),'ko','MarkerFaceColor','r')
    else
        plot(real(eig_A(i_pole)),imag(eig_A(i_pole)),'ko','MarkerFaceColor','g')
    end
    grid on
    xlabel('Re');
    ylabel('Im');
    subplot(2,1,2)
    bar(abs(participation_factor(i_pole,:)))
    xlim([0 pole_number+1])
	set(gca,'xTick',1:pole_number);
    set(gca,'XTickLabel',{'ild' 'ilq' 'vod' 'voq' 'iod' 'ioq' 'Gammad' 'Gammaq' 'Phid' 'Phiq' 'GammaSd' 'GammaSq' 'vplld' 'vpllq' 'Epsilonpll' 'DeltaThetapll' })
    grid on
    ylabel('Participation Factors');
    % end
end