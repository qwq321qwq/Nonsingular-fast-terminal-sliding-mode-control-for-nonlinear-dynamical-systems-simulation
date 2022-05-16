clear; clc; close all;


Kp = 1; Kd = 0.001; Kk = 0.01;


ddqa = [-0. -0.]'; dqa = [-0.1 -0.1]'; qa = [-0. 0.]';


g=9.81;
r1=1;r2=0.8;J1=5;J2=5;m1=0.5;m2=1.5;

q1=[];q2=[];dq1=[];dq2=[];ddq1=[];ddq2=[];
q=[];dq=[0;0];ddq=[];
qd1=[];qd2=[];dqd1=[];dqd2=[];ddqd1=[];ddqd2=[];
q1(1)=3.0;q2(1)=2.5;dq1(1)=0; dq2(1)=0; q(1,1)=3.0;q(2,1)=2.5;

for t = 1:700
    
    qd1(t)  =1.25-7/5*exp(-t)+7/20*exp(-4*t);qd2(t) =1.25+exp(-t)-(1/4)*exp(-4*t);
    qd      =[qd1(t);qd2(t)]; 
    dqd1(t) =7/5*exp(-t)-7/5*exp(-4*t);dqd2(t) =-exp(-t)+exp(-4*t);
    dqd     =[dqd1(t);dqd2(t)];
    ddqd1(t)=-7/5*exp(-t)+28/5*exp(-4*t);ddqd2(t) =exp(-t)-4*exp(-4*t);
    ddqd    =[ddqd1(t);ddqd2(t)];
    
    tau_d(:,t) =[2*sin(t)+0.5*sin(200*t);cos(2*t)+0.5*sin(200*t)];
    e11(t) = q1-qd1(t); e12(t) = q2-qd2(t); e1=[e11(t);e12(t)];
    e21(t) = dq1-dqd1(t); e22(t) = dq2-dqd2(t); e2=[e21(t);e22(t)];
    
    D=[(m1+m2)*r1^2+m2*r2^2+2*m2*r1*r2*cos(q2)+J1 m2*r2^2+m2*r1*r2*cos(q2);
        m2*r2^2+m2*r1*r2*cos(q2) m2*r2^2+J2];
    C=[-m2*r2^2+m2*r1*r2*sin(q2)*dq2^2-2*m2*r2^2+m2*r1*r2*sin(q2)*dq1*dq2;
       m2*r2^2+m2*r1*r2*sin(q2)*dq1^2];
    G=[(m1+m2)*g*r1*cos(q1)+m2*g*r2*cos(q1+q2);
        m2*g*r2*cos(q1+q2)];
        
    F2=-inv(D)*(C+G)-ddqd;
    
    E1=[1 0; 0 1];b0=12;b1=2.2;b2=2.8;
    gamma1=2; Gamma1=[gamma1 0;0 gamma1]; 
    gamma2=5/3; Gamma2=[gamma2 0; 0 gamma2];
    M1=1;M2=2;
%   sliding mode control
    rou = norm(inv(D))*(b0+b1*norm(q)+b2*norm(dq)^2);
    S = e1+[signa(e11(t),gamma1);signa(e12(t),gamma1)]+...
        [signa(e21(t),gamma2);signa(e22(t),gamma2)];
    
    Tau = -D*(M2*S+(rou+M1)*S/norm(S)+F2+inv(Gamma2)*...
        (E1+Gamma1*diag([signa(e11(t),(gamma1-1));signa(e12(t),(gamma1-1))]))*...
        [signa(e21(t),(2-gamma2));signa(e22(t),(2-gamma2))]);
    
    
%     所以下边反推了个寂寞
    deata_t = 0.01;
    ddq = inv(D) * (Tau - C - G + tau_d(:,t));
    ddq1= ddq(1); ddq2= ddq(2);
    dq  = ddq * deata_t + dq;
    dq1 = dq(1); dq2= dq(2);
    q   = dq * deata_t + q;
    q1  = q(1); q2= q(2);
    Q(:,t)=q;
    tauI(t,1)= Tau(1);
    tauI(t,2)= Tau(2);
%     Qd(qd1(t)

end

figure(1);
subplot(211);
plot(Q(1,:));
hold on;
plot(qd1);
% ylim([-0.3,0]);
% xlim([-0,1.5]);
xlabel('Test sample in joint1','FontSize',10);
% ylabel('Current','FontSize',10);
legend({'Tracing with controller','joint1'},'FontSize',10);

subplot(212);
plot(Q(2,:),'linewidth',1);
hold on;
plot(qd2,'linewidth',1);
% ylim([-0.3,0]);
xlabel('Test sample in joint2','FontSize',10);
% ylabel('Current','FontSize',10);
legend({'Tracing with controller','joint2'},'FontSize',10);


figure(2);
subplot(211);
plot(tauI(:,1));
hold on;
plot(tau_d(1,:)');
% ylim([-0.3,0]);
% xlim([-0,1.5]);
xlabel('Test sample in joint1','FontSize',10);
ylabel('Current','FontSize',10);
legend({'Tracing with controller','Torque1'},'FontSize',10);

subplot(212);
plot(tauI(:,2),'linewidth',1);
hold on;
plot(tau_d(2,:),'linewidth',1);
% ylim([-0.3,0]);
xlabel('Test sample in joint2','FontSize',10);
ylabel('Current','FontSize',10);
legend({'Tracing with controller','Torque2'},'FontSize',10);


function output = signa(x,a)
output =sign(x)*(norm(x,1)^a);
% output =diag(sign(x))*(norm(x,1)^a);
% output =(1./(1+exp(-x)))^v;
end
