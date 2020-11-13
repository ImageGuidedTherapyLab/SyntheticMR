

figure; hold on;
plot(TDpT1,squeeze(predstats(1,1,:)));
plot(TDpT1,squeeze(predstats(1,2,:)));
plot(TDpT1,squeeze(predstats(1,3,:)));
legend('GM','WM','CSF');
title('Mean Predicted T1 Value'); xlabel('Acquisition Spacing (s)'); ylabel('Mean T1 (s)');

figure; hold on;
plot(TDpT1,squeeze(predstats(2,1,:)));
plot(TDpT1,squeeze(predstats(2,2,:)));
plot(TDpT1,squeeze(predstats(2,3,:)));
legend('GM','WM','CSF');
title('Standard Deviation of Predicted T1 Value'); xlabel('Acquisition Spacing (s)'); ylabel('Std T1 (s)');

figure; plotyy(TDpT1,squeeze(predstats(2,1,:)),TDpT1,squeeze(predstats(2,2,:)));
title('Standard Deviation of Predicted T1 Value'); xlabel('Acquisition Spacing (s)'); ylabel('Std T1 (s)');
legend('GM','WM');

figure; hold on;
plot(TDpT1,squeeze(predstats(3,1,:)));
plot(TDpT1,squeeze(predstats(3,2,:)));
plot(TDpT1,squeeze(predstats(3,3,:)));
legend('GM','WM','CSF');
title('Mean Predicted T2 Value'); xlabel('Acquisition Spacing (s)'); ylabel('Mean T2 (s)');

figure; hold on;
plot(TDpT1,squeeze(predstats(4,1,:)));
plot(TDpT1,squeeze(predstats(4,2,:)));
plot(TDpT1,squeeze(predstats(4,3,:)));
legend('GM','WM','CSF');
title('Standard Deviation of Predicted T2 Value'); xlabel('Acquisition Spacing (s)'); ylabel('Std T2 (s)');

figure; plotyy(TDpT1,squeeze(predstats(4,1,:)),TDpT1,squeeze(predstats(4,2,:)));
title('Standard Deviation of Predicted T2 Value'); xlabel('Acquisition Spacing (s)'); ylabel('Std T2 (s)');
legend('GM','WM');

figure; hold on;
plot(TDpT1,squeeze(predstats(5,1,:)));
plot(TDpT1,squeeze(predstats(5,2,:)));
plot(TDpT1,squeeze(predstats(5,3,:)));
legend('GM','WM','CSF');
title('Mean Predicted M0 Value'); xlabel('Acquisition Spacing (s)'); ylabel('Mean M0 (s)');

figure; hold on;
plot(TDpT1,squeeze(predstats(6,1,:)));
plot(TDpT1,squeeze(predstats(6,2,:)));
plot(TDpT1,squeeze(predstats(6,3,:)));
legend('GM','WM','CSF');
title('Standard Deviation of Predicted M0 Value'); xlabel('Acquisition Spacing (s)'); ylabel('Std M0 (s)');
