function [] = energy_ratio(Y,mu,mass)

global n 
n = 15;

pos = Y(:,1:3*n);
vel = Y(:,3*n+1:6*n);

apo_pos = Y(:,1:3);

kinetic = zeros(length(apo_pos),n);
potential = zeros(length(apo_pos),n);
total = zeros(length(apo_pos),n);
ratio = zeros(length(apo_pos),n);

mu_sun = mu(2,1);

for i = 1:n
    for j = 1:length(apo_pos)
    if i == 2
        kinetic(:,i) = 0;
        potential(:,i) = 0;
        total(:,i) = kinetic(:,i) + potential(:,i); 
    else
        mass_i = mass(i,1);
        kinetic(j,i) = 1/2*mass_i*norm(vel(j,3*(i-1)+1:3*(i-1)+3))^2;
        potential(j,i) = -mass_i*mu_sun/norm(pos(j,3*(i-1)+1:3*(i-1)+3));
        total(j,i) = kinetic(j,i) + potential(j,i); 
        ratio(j,i) = kinetic(j,i)/potential(j,i);
    end
    end
end
time = linspace(1,length(pos),length(pos));

TOTAL = sum(total');
KINETIC = sum(kinetic');
POTENTIAL = sum(potential');
RATIO = KINETIC./POTENTIAL;

% l = 8;
% figure;
% subplot(l,2,[1 2])
% set(gca,'FontSize',14)
% plot(time,TOTAL)
% set(gca,'FontSize',14)
% title('Total energy: Solar System')
% axis([time(1) time(end) min(TOTAL) max(TOTAL)])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,3)
% plot(time,total(:,1))
% YTick_apo = [min(total(:,1));max(total(:,1))];
% set(gca,'FontSize',14,'YTick',YTick_apo)
% title('Apophis')
% axis([time(1) time(end) min(total(:,1)) max(total(:,1))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,4)
% plot(time,total(:,3))
% YTick_mer = [min(total(:,3));max(total(:,3))];
% set(gca,'FontSize',14,'YTick',YTick_mer)
% title('Mercury')
% axis([time(1) time(end) min(total(:,3)) max(total(:,3))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,5)
% plot(time,total(:,4))
% YTick_ven = [min(total(:,4));max(total(:,4))];
% set(gca,'FontSize',14,'YTick',YTick_ven)
% title('Venus')
% axis([time(1) time(end) min(total(:,4)) max(total(:,4))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,6)
% plot(time,total(:,5))
% YTick_ear = [min(total(:,5));max(total(:,5))];
% set(gca,'FontSize',14,'YTick',YTick_ear)
% title('Earth')
% axis([time(1) time(end) min(total(:,5)) max(total(:,5))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,7)
% plot(time,total(:,6))
% YTick_mar = [min(total(:,6));max(total(:,6))];
% set(gca,'FontSize',14,'YTick',YTick_mar)
% title('Mars')
% axis([time(1) time(end) min(total(:,6)) max(total(:,6))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,8)
% plot(time,total(:,7))
% YTick_jup = [min(total(:,7));max(total(:,7))];
% set(gca,'FontSize',14,'YTick',YTick_jup)
% title('Jupiter')
% axis([time(1) time(end) min(total(:,7)) max(total(:,7))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,9)
% plot(time,total(:,8))
% YTick_sat = [min(total(:,8));max(total(:,8))];
% set(gca,'FontSize',14,'YTick',YTick_sat)
% title('Saturn')
% axis([time(1) time(end) min(total(:,8)) max(total(:,8))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,10)
% plot(time,total(:,9))
% YTick_ura = [min(total(:,9));max(total(:,9))];
% set(gca,'FontSize',14,'YTick',YTick_ura)
% title('Uranus')
% axis([time(1) time(end) min(total(:,9)) max(total(:,9))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,11)
% plot(time,total(:,10))
% YTick_nep = [min(total(:,10));max(total(:,10))];
% set(gca,'FontSize',14,'YTick',YTick_nep)
% title('Neptune')
% axis([time(1) time(end) min(total(:,10)) max(total(:,10))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,12)
% plot(time,total(:,11))
% YTick_plu = [min(total(:,11));max(total(:,11))];
% set(gca,'FontSize',14,'YTick',YTick_plu)
% title('Pluto')
% axis([time(1) time(end) min(total(:,11)) max(total(:,11))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,13)
% plot(time,total(:,12))
% YTick_moo = [min(total(:,12));max(total(:,12))];
% set(gca,'FontSize',14,'YTick',YTick_moo)
% title('Moon')
% axis([time(1) time(end) min(total(:,12)) max(total(:,12))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,14)
% plot(time,total(:,13))
% YTick_cer = [min(total(:,13));max(total(:,13))];
% set(gca,'FontSize',14,'YTick',YTick_cer)
% title('Ceres')
% axis([time(1) time(end) min(total(:,13)) max(total(:,13))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,15)
% plot(time,total(:,14))
% YTick_pal = [min(total(:,14));max(total(:,14))];
% set(gca,'FontSize',14,'YTick',YTick_pal)
% title('Pallas')
% axis([time(1) time(end) min(total(:,14)) max(total(:,14))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,16)
% plot(time,total(:,15))
% YTick_ves = [min(total(:,15));max(total(:,15))];
% set(gca,'FontSize',14,'YTick',YTick_ves)
% title('Vesta')
% axis([time(1) time(end) min(total(:,15)) max(total(:,15))])
% set(gca,'xtick',[], 'xticklabel',{})



% l = 8;
% figure;
% subplot(l,2,[1 2])
% set(gca,'FontSize',14)
% plot(time,RATIO)
% title('Ratio of kinetic over potential: Solar System')
% axis([time(1) time(end) min(RATIO) max(RATIO)])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,3)
% plot(time,ratio(:,1))
% YTick_apo = [min(ratio(:,1));max(ratio(:,1))];
% set(gca,'FontSize',14,'YTick',YTick_apo)
% title('Apophis')
% axis([time(1) time(end) min(ratio(:,1)) max(ratio(:,1))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,4)
% plot(time,ratio(:,3))
% YTick_mer = [min(ratio(:,3));max(ratio(:,3))];
% set(gca,'FontSize',14,'YTick',YTick_mer)
% title('Mercury')
% axis([time(1) time(end) min(ratio(:,3)) max(ratio(:,3))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,5)
% plot(time,ratio(:,4))
% YTick_ven = [min(ratio(:,4));max(ratio(:,4))];
% set(gca,'FontSize',14,'YTick',YTick_ven)
% title('Venus')
% axis([time(1) time(end) min(ratio(:,4)) max(ratio(:,4))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,6)
% plot(time,ratio(:,5))
% YTick_ear = [min(ratio(:,5));max(ratio(:,5))];
% set(gca,'FontSize',14,'YTick',YTick_ear)
% title('Earth')
% axis([time(1) time(end) min(ratio(:,5)) max(ratio(:,5))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,7)
% plot(time,ratio(:,6))
% YTick_mar = [min(ratio(:,6));max(ratio(:,6))];
% set(gca,'FontSize',14,'YTick',YTick_mar)
% title('Mars')
% axis([time(1) time(end) min(ratio(:,6)) max(ratio(:,6))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,8)
% plot(time,ratio(:,7))
% YTick_jup = [min(ratio(:,7));max(ratio(:,7))];
% set(gca,'FontSize',14,'YTick',YTick_jup)
% title('Jupiter')
% axis([time(1) time(end) min(ratio(:,7)) max(ratio(:,7))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,9)
% plot(time,ratio(:,8))
% YTick_sat = [min(ratio(:,8));max(ratio(:,8))];
% set(gca,'FontSize',14,'YTick',YTick_sat)
% title('Saturn')
% axis([time(1) time(end) min(ratio(:,8)) max(ratio(:,8))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,10)
% plot(time,ratio(:,9))
% YTick_ura = [min(ratio(:,9));max(ratio(:,9))];
% set(gca,'FontSize',14,'YTick',YTick_ura)
% title('Uranus')
% axis([time(1) time(end) min(ratio(:,9)) max(ratio(:,9))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,11)
% plot(time,ratio(:,10))
% YTick_nep = [min(ratio(:,10));max(ratio(:,10))];
% set(gca,'FontSize',14,'YTick',YTick_nep)
% title('Neptune')
% axis([time(1) time(end) min(ratio(:,10)) max(ratio(:,10))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,12)
% plot(time,ratio(:,11))
% YTick_plu = [min(ratio(:,11));max(ratio(:,11))];
% set(gca,'FontSize',14,'YTick',YTick_plu)
% title('Pluto')
% axis([time(1) time(end) min(ratio(:,11)) max(ratio(:,11))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,13)
% plot(time,ratio(:,12))
% YTick_moo = [min(ratio(:,12));max(ratio(:,12))];
% set(gca,'FontSize',14,'YTick',YTick_moo)
% title('Moon')
% axis([time(1) time(end) min(ratio(:,12)) max(ratio(:,12))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,14)
% plot(time,ratio(:,13))
% YTick_cer = [min(ratio(:,13));max(ratio(:,13))];
% set(gca,'FontSize',14,'YTick',YTick_cer)
% title('Ceres')
% axis([time(1) time(end) min(ratio(:,13)) max(ratio(:,13))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,15)
% plot(time,ratio(:,14))
% YTick_pal = [min(ratio(:,14));max(ratio(:,14))];
% set(gca,'FontSize',14,'YTick',YTick_pal)
% set(gca,'FontSize',14)
% title('Pallas')
% axis([time(1) time(end) min(ratio(:,14)) max(ratio(:,14))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,16)
% plot(time,ratio(:,15))
% YTick_ves = [min(ratio(:,15));max(ratio(:,15))];
% set(gca,'FontSize',14,'YTick',YTick_ves)
% title('Vesta')
% axis([time(1) time(end) min(ratio(:,15)) max(ratio(:,15))])
% set(gca,'xtick',[], 'xticklabel',{})



mean_total(:,1) = mean(total(:,1));
mean_total(:,2) = mean(total(:,2));
mean_total(:,3) = mean(total(:,3));
mean_total(:,4) = mean(total(:,4));
mean_total(:,5) = mean(total(:,5));
mean_total(:,6) = mean(total(:,6));
mean_total(:,7) = mean(total(:,7));
mean_total(:,8) = mean(total(:,8));
mean_total(:,9) = mean(total(:,9));
mean_total(:,10) = mean(total(:,10));
mean_total(:,11) = mean(total(:,11));
mean_total(:,12) = mean(total(:,12));
mean_total(:,13) = mean(total(:,13));
mean_total(:,14) = mean(total(:,14));
mean_total(:,15) = mean(total(:,15));
mean_TOTAL = mean(TOTAL);

% mean_TOTAL = mean(TOTAL);
% % mean_kinetic = mean(kinetic);
% % mean_potential = mean(potential);
% total = total/mean_total;
% TOTAL = TOTAL/mean_TOTAL;

for j = 1:15
    total(:,j) = total(:,j)/mean_total(:,j);
end

TOTAL = TOTAL/mean_TOTAL;
% 
% l = 8;
% figure;
% subplot(l,2,[1 2])
% set(gca,'FontSize',14)
% plot(time,TOTAL)
% title('Normalized total energy: Solar System')
% axis([time(1) time(end) min(TOTAL) max(TOTAL)])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,3)
% plot(time,total(:,1))
% set(gca,'FontSize',14)
% title('Apophis')
% axis([time(1) time(end) min(total(:,1)) max(total(:,1))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,4)
% plot(time,total(:,3))
% set(gca,'FontSize',14)
% title('Mercury')
% axis([time(1) time(end) min(total(:,3)) max(total(:,3))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,5)
% plot(time,total(:,4))
% set(gca,'FontSize',14)
% title('Venus')
% axis([time(1) time(end) min(total(:,4)) max(total(:,4))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,6)
% plot(time,total(:,5))
% set(gca,'FontSize',14)
% title('Earth')
% axis([time(1) time(end) min(total(:,5)) max(total(:,5))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,7)
% plot(time,total(:,6))
% set(gca,'FontSize',14)
% title('Mars')
% axis([time(1) time(end) min(total(:,6)) max(total(:,6))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,8)
% plot(time,total(:,7))
% set(gca,'FontSize',14)
% title('Jupiter')
% axis([time(1) time(end) min(total(:,7)) max(total(:,7))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,9)
% plot(time,total(:,8))
% set(gca,'FontSize',14)
% title('Saturn')
% axis([time(1) time(end) min(total(:,8)) max(total(:,8))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,10)
% plot(time,total(:,9))
% set(gca,'FontSize',14)
% title('Uranus')
% axis([time(1) time(end) min(total(:,9)) max(total(:,9))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,11)
% plot(time,total(:,10))
% set(gca,'FontSize',14)
% title('Neptune')
% axis([time(1) time(end) min(total(:,10)) max(total(:,10))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,12)
% plot(time,total(:,11))
% set(gca,'FontSize',14)
% title('Pluto')
% axis([time(1) time(end) min(total(:,11)) max(total(:,11))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,13)
% plot(time,total(:,12))
% set(gca,'FontSize',14)
% title('Moon')
% axis([time(1) time(end) min(total(:,12)) max(total(:,12))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,14)
% plot(time,total(:,13))
% set(gca,'FontSize',14)
% title('Ceres')
% axis([time(1) time(end) min(total(:,13)) max(total(:,13))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,15)
% plot(time,total(:,14))
% set(gca,'FontSize',14)
% title('Pallas')
% axis([time(1) time(end) min(total(:,14)) max(total(:,14))])
% set(gca,'xtick',[], 'xticklabel',{})
% subplot(l,2,16)
% plot(time,total(:,15))
% set(gca,'FontSize',14)
% title('Vesta')
% axis([time(1) time(end) min(total(:,15)) max(total(:,15))])
% set(gca,'xtick',[], 'xticklabel',{})

l = 8;
%figure;
subplot(l,2,[1 2])
set(gca,'FontSize',14)
plot(time,TOTAL)
set(gca,'FontSize',14)
title('Normalized total energy: Solar System')
axis([time(1) time(end) min(TOTAL) max(TOTAL)])
set(gca,'xtick',[], 'xticklabel',{})
subplot(l,2,3)
plot(time,total(:,1))
YTick_apo = [min(total(:,1));max(total(:,1))];
set(gca,'FontSize',14,'YTick',YTick_apo)
title('Apophis')
axis([time(1) time(end) min(total(:,1)) max(total(:,1))])
set(gca,'xtick',[], 'xticklabel',{})
subplot(l,2,4)
plot(time,total(:,3))
YTick_mer = [min(total(:,3));max(total(:,3))];
set(gca,'FontSize',14,'YTick',YTick_mer)
title('Mercury')
axis([time(1) time(end) min(total(:,3)) max(total(:,3))])
set(gca,'xtick',[], 'xticklabel',{})
subplot(l,2,5)
plot(time,total(:,4))
YTick_ven = [min(total(:,4));max(total(:,4))];
set(gca,'FontSize',14,'YTick',YTick_ven)
title('Venus')
axis([time(1) time(end) min(total(:,4)) max(total(:,4))])
set(gca,'xtick',[], 'xticklabel',{})
subplot(l,2,6)
plot(time,total(:,5))
YTick_ear = [min(total(:,5));max(total(:,5))];
set(gca,'FontSize',14,'YTick',YTick_ear)
title('Earth')
axis([time(1) time(end) min(total(:,5)) max(total(:,5))])
set(gca,'xtick',[], 'xticklabel',{})
subplot(l,2,7)
plot(time,total(:,6))
YTick_mar = [min(total(:,6));max(total(:,6))];
set(gca,'FontSize',14,'YTick',YTick_mar)
title('Mars')
axis([time(1) time(end) min(total(:,6)) max(total(:,6))])
set(gca,'xtick',[], 'xticklabel',{})
subplot(l,2,8)
plot(time,total(:,7))
YTick_jup = [min(total(:,7));max(total(:,7))];
set(gca,'FontSize',14,'YTick',YTick_jup)
title('Jupiter')
axis([time(1) time(end) min(total(:,7)) max(total(:,7))])
set(gca,'xtick',[], 'xticklabel',{})
subplot(l,2,9)
plot(time,total(:,8))
YTick_sat = [min(total(:,8));max(total(:,8))];
set(gca,'FontSize',14,'YTick',YTick_sat)
title('Saturn')
axis([time(1) time(end) min(total(:,8)) max(total(:,8))])
set(gca,'xtick',[], 'xticklabel',{})
subplot(l,2,10)
plot(time,total(:,9))
YTick_ura = [min(total(:,9));max(total(:,9))];
set(gca,'FontSize',14,'YTick',YTick_ura)
title('Uranus')
axis([time(1) time(end) min(total(:,9)) max(total(:,9))])
set(gca,'xtick',[], 'xticklabel',{})
subplot(l,2,11)
plot(time,total(:,10))
YTick_nep = [min(total(:,10));max(total(:,10))];
set(gca,'FontSize',14,'YTick',YTick_nep)
title('Neptune')
axis([time(1) time(end) min(total(:,10)) max(total(:,10))])
set(gca,'xtick',[], 'xticklabel',{})
subplot(l,2,12)
plot(time,total(:,11))
YTick_plu = [min(total(:,11));max(total(:,11))];
set(gca,'FontSize',14,'YTick',YTick_plu)
title('Pluto')
axis([time(1) time(end) min(total(:,11)) max(total(:,11))])
set(gca,'xtick',[], 'xticklabel',{})
subplot(l,2,13)
plot(time,total(:,12))
YTick_moo = [min(total(:,12));max(total(:,12))];
set(gca,'FontSize',14,'YTick',YTick_moo)
title('Moon')
axis([time(1) time(end) min(total(:,12)) max(total(:,12))])
set(gca,'xtick',[], 'xticklabel',{})
subplot(l,2,14)
plot(time,total(:,13))
YTick_cer = [min(total(:,13));max(total(:,13))];
set(gca,'FontSize',14,'YTick',YTick_cer)
title('Ceres')
axis([time(1) time(end) min(total(:,13)) max(total(:,13))])
set(gca,'xtick',[], 'xticklabel',{})
subplot(l,2,15)
plot(time,total(:,14))
YTick_pal = [min(total(:,14));max(total(:,14))];
set(gca,'FontSize',14,'YTick',YTick_pal)
title('Pallas')
axis([time(1) time(end) min(total(:,14)) max(total(:,14))])
set(gca,'xtick',[], 'xticklabel',{})
subplot(l,2,16)
plot(time,total(:,15))
YTick_ves = [min(total(:,15));max(total(:,15))];
set(gca,'FontSize',14,'YTick',YTick_ves)
title('Vesta')
axis([time(1) time(end) min(total(:,15)) max(total(:,15))])
set(gca,'xtick',[], 'xticklabel',{})


end




