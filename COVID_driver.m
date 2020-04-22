V0 = 0.02;
H0 = 1;
L0 = 0;
I0 = 0;
M0 = 0;
F0 = 0;
R0 = 0;
E0 = 1;
P0 = 1;
A0 = 1;
S0 = 0.1;

x0 = [V0 H0 L0 I0 M0 F0 R0 E0 P0 A0 S0];
pars = [];

[t,x] = ode15s(@(t,x) hancioglu(t,x,pars), [0 25], x0);
size(x)
% extract variables
V = x(:,1)*1.5e6;
H = x(:,2);
L = x(:,3);
I = x(:,4);
M = x(:,5);
F = x(:,6);
R = x(:,7);
E = x(:,8);
P = x(:,9);
A = x(:,10);
S = x(:,11);
D = 1-H-R-I-L;

figure(1), clf
subplot(4,3,1)
semilogy(t,V)
ylabel('V');
subplot(4,3,2)
plot(t,H)
ylabel('H');
subplot(4,3,3)
plot(t,I)
ylabel('I');
subplot(4,3,4)
plot(t,M)
ylabel('M');
subplot(4,3,5)
semilogy(t,F)
ylabel('F');
subplot(4,3,6)
plot(t,R)
ylabel('R');
subplot(4,3,7)
semilogy(t,E)
ylabel('E');
subplot(4,3,8)
semilogy(t,P)
ylabel('P');
subplot(4,3,9)
semilogy(t,A)
ylabel('A');
subplot(4,3,10)
plot(t,D)
ylabel('D');
subplot(4,3,11)
plot(t,S)
ylabel('S');
subplot(4,3,12)
plot(t,L)
ylabel('L');

% FIGURE 2
wolfel_sputum = load('wolfel_sputum2.dat');
yan_throat    = load('yan_throat.dat');
yan_sputum    = load('yan_sputum.dat');

figure(2), clf
subplot(1,2,1)
semilogy(t-7,V,'-k','Linewidth',2)
hold on
scatter(yan_sputum(:,1),yan_sputum(:,2)*5e3,'ob','Linewidth',4)
scatter(yan_throat(:,1),yan_throat(:,2)*5e3,'sr','Linewidth',4)
set(gca,'Fontsize',18);
xlim([0 10])
ylim([10^2 10^9])
ylabel('Viral load /mL','Fontsize',22)
xlabel('Days','Fontsize',22)
legend({'Model','Sputum','Throat'},'Fontsize',18,'Location','southwest')
text(0.5,5e8,'A','Fontsize',17)

subplot(1,2,2)
semilogy(t-7,V,'-k','Linewidth',2)
hold on
scatter(wolfel_sputum(:,1),wolfel_sputum(:,2),'sg','Linewidth',4)
set(gca,'Fontsize',16);
xlim([0 10])
ylim([10^2 10^9])
ylabel('Viral load /mL','Fontsize',22)
xlabel('Days','Fontsize',22)
text(0.5,5e8,'B','Fontsize',18)

% FIGURE 3
figure(3), clf
subplot(3,2,1)
plot(t,H,'Linewidth',2)
hold on
plot(t,L,'Linewidth',2)
plot(t,I,'Linewidth',2)
plot(t,R,'Linewidth',2)
plot(t,D,'Linewidth',2)
set(gca,'Fontsize',16);
legend({'Healthy (H)','Latent (L)','Infected (I)','Resistant (R)','Damaged (D)'},'Location','southwest','Fontsize',15)
ylabel('Fractional population','Fontsize',18)
xlabel('Days','Fontsize',18)
ylim([0 1])
xlim([0 25])
text(0.5,0.95,'A','Fontsize',17)

subplot(3,2,2)
semilogy(t,F,'-k','Linewidth',2)
set(gca,'Fontsize',16);
ylabel('Interferon (F)','Fontsize',18)
xlabel('Days','Fontsize',18)
xlim([0 25])
text(0.5,2e4,'B','Fontsize',17)

subplot(3,2,3)
plot(t, M,'-k','Linewidth',2)
set(gca,'Fontsize',16);
ylabel('Antigen presenting cells (M)','Fontsize',18)
xlabel('Days','Fontsize',18)
ylim([0 1])
xlim([0 25])
text(0.5,0.95,'C','Fontsize',17)

subplot(3,2,4)
semilogy(t, E,'-k','Linewidth',2)
set(gca,'Fontsize',16);
ylabel('Effector cells (E)','Fontsize',18)
xlabel('Days','Fontsize',18)
ylim([0 10^2])
xlim([0 25])
text(0.5,0.8e2,'D','Fontsize',17)

subplot(3,2,5)
semilogy(t, A,'-k','Linewidth',2)
set(gca,'Fontsize',16);
ylabel('Antibodies (A)','Fontsize',18)
xlabel('Days','Fontsize',18)
ylim([10^-5 10^4])
xlim([0 25])
yticks([1e-4 1e-2 1 1e2 1e4])
text(0.5,3e3,'E','Fontsize',17)

subplot(3,2,6)
plot(t, S,'-k','Linewidth',2)
set(gca,'Fontsize',16);
ylabel('Specificity (S)','Fontsize',18)
xlabel('Days','Fontsize',18)
ylim([0 1])
xlim([0 25])
text(0.5,0.95,'F','Fontsize',17)
