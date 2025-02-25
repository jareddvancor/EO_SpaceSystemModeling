sun = readmatrix('SolarSpec.txt')
qe = [[350:50:1050]',[0.3,0.6,0.8,0.95,0.96,0.94,0.93,0.92,0.85,0.75,0.6,0.45,0.27,0.1,0.05]']
opt = [qe(:,1),qe(:,1)*0+.95]
tgt = [qe(:,1),qe(:,1)*0+.3]

plot(qe(:,1),qe(:,2)); hold on
sun2 = [[350:1050]', interp1(sun(:,1),sun(:,2),[350:1050])']
plot(sun2(:,1),sun2(:,2)/max(sun2(:,2)))
plot(opt(:,1),opt(:,2),'-');
plot(tgt(:,1),tgt(:,2),'-');
xline(400,'--');
xline(1000,'--')
legend({'Quantum Efficency','Normalize Sun','Optical Transmission','Target Reflectance','cuton','cutoff'})
title('Relative Spectral Sensitivity')
xlabel('Wavelength [nm]')
ylabel('Relative Sensitivivity []')


sun3 = interp1(sun2(:,1),sun2(:,2),[400:1000]);
qe2 = interp1(qe(:,1),qe(:,2),[400:1000])
cwl = mean([400:1000].*sun3/mean(sun3).*qe2/mean(qe2))


