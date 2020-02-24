alphaKQ = data_alpha; thetaKQ = data_theta; VmKQ = data_Vm;


%%
alphaKinf = data_alpha; thetaKinf = data_theta; VmKinf = data_Vm;

t = alphaKQ(:,1);
ref = thetaKQ(:,2);

figure(1)
clf
subplot(3,2,1) % Alpha
plot(t,alphaKQ(:,2),  'b')
hold on
plot(t,alphaKinf(:,2),'r')
title('Alpha (\alpha)','FontSize', 16)
legend('\alpha_{K_Q}','\alpha_{K_{H\infty}}','Interpreter','latex')
axis([0 20 -15 15])
set(gca,'XGrid','on', 'YGrid', 'on', 'FontSize', 16)

subplot(3,2,2) % Alpha
plot(t,alphaKQ(:,2),  'b')
hold on
plot(t,alphaKinf(:,2),'r')
title('Alpha (\alpha)','FontSize', 16)
legend('\alpha_{K_Q}','\alpha_{K_{H\infty}}','Interpreter','latex')
axis([9.5 13 -15 15])
set(gca,'XGrid','on', 'YGrid', 'on', 'FontSize', 16)

subplot(3,2,3) % Theta
plot(t,thetaKQ(:,3),  'b')
hold on
plot(t,thetaKinf(:,3),'r')
plot(t,ref,'k') % Referência
title('Theta (\theta)','FontSize', 16)
legend('\theta_{K_Q}','\theta_{K_{H\infty}}','\theta_{ref}','Interpreter','latex')
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 16)

subplot(3,2,4) % Theta
plot(t,thetaKQ(:,3),  'b')
hold on
plot(t,thetaKinf(:,3),'r')
plot(t,ref,'k') % Referência
title('Theta (\theta)','FontSize', 16)
legend('\theta_{K_Q}','\theta_{K_{H\infty}}','\theta_{ref}','Interpreter','latex')
axis([9.5 13 -20 20])
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 16)

subplot(3,2,5) % Vm
plot(t,VmKQ(:,2),  'b')
hold on
plot(t,VmKinf(:,2),'r')
title('Vm','FontSize', 16)
legend('Vm_{K_Q}','Vm_{K_{H\infty}}','Interpreter','latex')
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 16)
axis([0 20 -5 max(VmKinf(:,2))])

subplot(3,2,6) % Vm
plot(t,VmKQ(:,2),  'b')
hold on
plot(t,VmKinf(:,2),'r')
title('Vm','FontSize', 16)
legend('Vm_{K_Q}','Vm_{K_{H\infty}}','Interpreter','latex')
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 16)
axis([9.5 13 -5 max(VmKinf(:,2))])
%min(VmKinf(:,2))
%% Índices de desempenho
clc
MSEQ_theta = mse(ref(1200:end),thetaKQ((1200:end),3));
MSEinf_theta = mse(ref(1200:end),thetaKinf((1200:end),3));
disp('MSE theta- [KQ    Kinf]')
[MSEQ_theta MSEinf_theta]

RMSEQ_theta = sqrt(mean((ref(1200:end)-thetaKQ((1200:end),3)).^2));
RMSEinf_theta = sqrt(mean((ref(1200:end)-thetaKinf((1200:end),3)).^2));
disp('RMSE theta- [KQ    Kinf]')
[RMSEQ_theta RMSEinf_theta]
pause

ref_alpha = zeros(length(alphaKQ),1);
MSEQ_alpha = mse(ref_alpha(1200:end),alphaKQ((1200:end),2));
MSEinf_alpha = mse(ref_alpha(1200:end),alphaKinf((1200:end),2));
disp('MSE alpha - [KQ    Kinf]')
[MSEQ_alpha MSEinf_alpha]

RMSEQ_alpha = sqrt(mse(ref_alpha(1200:end),alphaKQ((1200:end),2)));
RMSEinf_alpha = sqrt(mse(ref_alpha(1200:end),alphaKinf((1200:end),2)));
disp('RMSE alpha - [KQ    Kinf]')
[RMSEQ_alpha RMSEinf_alpha]
pause


ref_Vm = ref_alpha;
MSEQ_Vm = mse(ref_Vm(1200:end),VmKQ((1200:end),2));
MSEinf_Vm = mse(ref_Vm(1200:end),VmKinf((1200:end),2));
disp('MSE Vm - [KQ    Kinf]')
[MSEQ_Vm MSEinf_Vm]

RMSEQ_Vm = sqrt(mse(ref_Vm(1200:end),VmKQ((1200:end),2)));
RMSEinf_Vm = sqrt(mse(ref_Vm(1200:end),VmKinf((1200:end),2)));
disp('RMSE Vm - [KQ    Kinf]')
[RMSEQ_Vm RMSEinf_Vm]
