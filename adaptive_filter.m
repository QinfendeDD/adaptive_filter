% %-----------LMS�㷨------------%
% close all;
% clear  all;
% clc;
% p=9;
% loop=5;
% N=1500;
% n=0:N-1;
% x=sin(600.1*pi*n);
% x1=0.8*sin(600.1*pi*n);
% s=1/6*sin(2000.1*pi*n);
% d=s+x1;
% 
% figure;
% subplot(5,1,1);
% plot(x);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('Զ�������ź�');
% 
% subplot(5,1,2);
% plot(x1);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% subplot(5,1,3);
% plot(s);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('���������ź�');
% 
% subplot(5,1,4);
% plot(d);  grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% y=zeros(size(n));
% w=zeros(1,p+1);
% e=zeros(size(n));
% E=zeros(size(n));
% u=0.001;
% for q=1:loop
%    for i=1:N-p
%     y(i+p)=w*x(i:i+p)';
%     e(i+p)=d(i+p)-y(i+p);
%     w=w+2*u*e(i+p)*x(i:i+p);
%    end
%    E=E+e.^2;
% end
%  E=E/loop;
%  subplot(5,1,5);
%  plot(y);  grid;
%  xlabel('����n');
%  ylabel('��ֵ');
%  title('����Ӧ�˲�����ź�');
% 
%  figure;
%  plot(E);  grid;
%  figure;
%  plot(d(p+1:N),'g');  grid;   
%  hold on;
%  plot(y(p+1:N),'b');   grid;
%  hold on;
%  ds=y-d;
%  plot(ds(p+1:N),'r');  grid;
%  legend('�����ź�','����Ӧ�˲�������ź�','���ֵ');



% %-------------�������Ƶı䲽��LMS�㷨----------%
% close all;
% clear  all;
% clc;
% p=9;
% loop=1;
% N=1500;
% a=0.97;    r=0.006; %�䲽���������Ʋ���
% n=0:N-1;
% x=sin(600.1*pi*n);
% x1=0.8*sin(600.1*pi*n);
% s=1/6*sin(2000.1*pi*n);
% d=s+x1;
% 
% figure;
% subplot(5,1,1);
% plot(x);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('Զ�������ź�');
% 
% subplot(5,1,2);
% plot(x1);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% subplot(5,1,3);
% plot(s);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('���������ź�');
% 
% subplot(5,1,4);
% plot(d);  grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% y=zeros(size(n));
% w=zeros(1,p+1);
% e=zeros(size(n));
% E=zeros(size(n));
% u=zeros(size(n));
% for q=1:loop
% for i=1:N-p
%     y(i+p)=w*x(i:i+p)';
%     e(i+p)=d(i+p)-y(i+p);
%     u(i+p)=a*u(i+p-1)+r*e(i+p-1)*e(i+p-1);
%     w=w+u(i+p)*e(i+p)*x(i:i+p);
% end
%   E=E+e.^2;
% end
%  E=E/loop;
%  subplot(5,1,5);
%  plot(y);  grid;
%  xlabel('����n');
%  ylabel('��ֵ');
%  title('����Ӧ�˲�����ź�');
% 
% 
%  figure;
%  plot(E);  grid;
%  
%  figure;
%  plot(u(p+1:N));  grid;
%  xlabel('����n');
%  ylabel('mu��ȡֵ');
%  title('�������ӱ仯����');
%  
%  figure;
%  plot(d(p+1:N),'g');  grid;   
%  hold on;
%  plot(y(p+1:N),'b');   grid;
%  hold on;
%  ds=y-d;
%  plot(ds(p+1:N),'r');  grid;
%  legend('�����ź�','����Ӧ�˲�������ź�','���ֵ');
% 




% %-------�����Ƶı䲽��LMS�㷨����ֹ������������ٶȹ���-------%
% close all;
% clear  all;
% clc;
% p=9;
% loop=1;
% N=500;
% a=0.97;    r=0.006;   %�䲽���������Ʋ���
% n=0:N-1;
% x=sin(600.1*pi*n);
% x1=0.8*sin(600.1*pi*n);
% s=1/6*sin(2000.1*pi*n);
% d=s+x1;
% 
% figure;
% subplot(5,1,1);
% plot(x);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('Զ�������ź�');
% 
% subplot(5,1,2);
% plot(x1);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% subplot(5,1,3);
% plot(s);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('���������ź�');
% 
% subplot(5,1,4);
% plot(d);  grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% y=zeros(size(n));
% w=zeros(1,p+1);
% e=zeros(size(n));
% E=zeros(size(n));
% u=zeros(size(n));
% u1=0.005;
% u2=0.00001;
% for q=1:loop
% for i=1:N-p
%     y(i+p)=w*x(i:i+p)';
%     e(i+p)=d(i+p)-y(i+p);
%     u(i+p)=a*u(i+p-1)+r*e(i+p)*e(i+p);
%     if u(i+p)>u1
%         u(i+p)=u1;
%      else if u(i+p)<u2
%         u(i+p)=u2;
%          end
%     end
%     w=w+u(i+p)*e(i+p)*x(i:i+p);
% end
%   E=E+e.^2;
% end
%  E=E/loop;
%  subplot(5,1,5);
%  plot(y);  grid;
%  xlabel('����n');
%  ylabel('��ֵ');
%  title('����Ӧ�˲�����ź�');
% 
% 
%  figure;
%  plot(E);  grid;
%  
%  figure;
%  plot(u(p+1:N));  grid;
%  xlabel('����n');
%  ylabel('mu��ȡֵ');
%  title('�������ӱ仯����');
%  
%  figure;
%  plot(d(p+1:N),'g');  grid;   
%  hold on;
%  plot(y(p+1:N),'b');   grid;
%  hold on;
%  ds=y-d;
%  plot(ds(p+1:N),'r');  grid;
%  legend('�����ź�','����Ӧ�˲�������ź�','���ֵ');



%----------�Ľ��ı䲽���������ƣ�LMS�㷨------------%
close all;
clear  all;
clc;
p=9;
loop=1;
N=1000;
a=0.97;    r=0.006;   %�䲽���������Ʋ���
n=0:N-1;
x=sin(600.1*pi*n);
x1=0.8*sin(600.1*pi*n);
s=1/6*sin(2000.1*pi*n);
d=s+x1;

figure;
subplot(5,1,1);
plot(x);   grid;
xlabel('����n');
ylabel('��ֵ');
title('Զ�������ź�');

subplot(5,1,2);
plot(x1);   grid;
xlabel('����n');
ylabel('��ֵ');
title('�����ź�');

subplot(5,1,3);
plot(s);   grid;
xlabel('����n');
ylabel('��ֵ');
title('���������ź�');

subplot(5,1,4);
plot(d);  grid;
xlabel('����n');
ylabel('��ֵ');
title('�����ź�');

y=zeros(size(n));
w=zeros(1,p+1);
e=zeros(size(n));
E=zeros(size(n));
u=zeros(size(n));
P=zeros(size(n));
u1=0.005;
u2=0.00001;
a=0.97;
r=0.006;
b=0.99;
for q=1:loop
for i=1:N-p
    y(i+p)=w*x(i:i+p)';
    e(i+p)=d(i+p)-y(i+p);
    P(i+p)=b*P(i+p-1)+(1-b)*e(i+p)*e(i+p-1);
    u(i+p)=a*u(i+p-1)+r*P(i+p);
    if u(i+p)>u1
        u(i+p)=u1;
     else if u(i+p)<u2
        u(i+p)=u2;
       end
    end
    w=w+u(i+p)*e(i+p)*x(i:i+p);
end
    E=E+e.^2;
end
 E=E/loop;
 subplot(5,1,5);
 plot(y);  grid;
 xlabel('����n');
 ylabel('��ֵ');
 title('����Ӧ�˲�����ź�');

 figure;
 plot(E);  grid;
 
 figure;
 plot(u(p+1:N));  grid;
 xlabel('����n');
 ylabel('mu��ȡֵ');
 title('�������ӱ仯����');
 
 figure;
 plot(d(p+1:N),'g');  grid;   
 hold on;
 plot(y(p+1:N),'b');   grid;
 hold on;
 ds=y-d;
 plot(ds(p+1:N),'r');  grid;
 legend('�����ź�','����Ӧ�˲�������ź�','���ֵ');


% 
% %------------NLMS-------------------%
% close all;
% clear  all;
% clc;
% p=9;
% loop=100;
% N=500;
% n=0:N-1;
% x=sin(600.1*pi*n);
% x1=0.8*sin(600.1*pi*n);
% s=1/6*sin(2000.1*pi*n);
% d=s+x1;
% 
% figure;
% subplot(5,1,1);
% plot(x);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('Զ�������ź�');
% 
% subplot(5,1,2);
% plot(x1);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% subplot(5,1,3);
% plot(s);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('���������ź�');
% 
% subplot(5,1,4);
% plot(d);  grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% y=zeros(size(n));
% w=zeros(1,p+1);
% e=zeros(size(n));
% E=zeros(size(n));
% u=1.0;%0<u<2
% b=1.0e-12;
% for q=1:loop
% for i=1:N-p
%     y(i+p)=w*x(i:i+p)';
%     e(i+p)=d(i+p)-y(i+p);
%     w=w+u/(b+x(i:i+p)*x(i:i+p).')*e(i+p)*x(i:i+p);
% end
% E=E+e.^2;
% end
%  E=E/loop;
%  subplot(5,1,5);
%  plot(y);  grid;
%  xlabel('����n');
%  ylabel('��ֵ');
%  title('����Ӧ�˲�����ź�');
% 
% 
%  figure;
%  plot(E(p+1:N));  grid;
%  
%  figure;
%  plot(d(p+1:N),'g');  grid;   
%  hold on;
%  plot(y(p+1:N),'b');   grid;
%  hold on;
%  ds=y-d;
%  plot(ds(p+1:N),'r');  grid;
%  legend('�����ź�','����Ӧ�˲�������ź�','���ֵ');
% 
% 
% %------------PNLMS---------------%
% close all;
% clear  all;
% clc;
% p=9;
% loop=100;
% N=500;
% n=0:N-1;
% x=sin(600.1*pi*n);
% x1=0.8*sin(600.1*pi*n);
% s=1/6*sin(2000.1*pi*n);
% d=s+x1;
% 
% figure;
% subplot(5,1,1);
% plot(x);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('Զ�������ź�');
% 
% subplot(5,1,2);
% plot(x1);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% subplot(5,1,3);
% plot(s);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('���������ź�');
% 
% subplot(5,1,4);
% plot(d);  grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% y=zeros(size(n));
% w=zeros(1,p+1);
% e=zeros(size(n));
% r=zeros(1,p+1);
% E=zeros(size(n));
% u=1.0;       sigma=1.0e-12;  
% ro=1/(p+1);  sigmap=0.01;
% for q=1:loop
% for i=1:N-p
%     y(i+p)=w*x(i:i+p)';
%     e(i+p)=d(i+p)-y(i+p);
%     temp=abs(w);
%     temp=max(temp);
%     rmax=ro*max(temp,sigmap);
%     for k=1:p+1
%         r(k)=max(rmax,w(k));
%     end
%     g=(p+1)*r./sum(abs(r));
%     g=diag(g);
%     w=w+u/(x(i:i+p)*g*x(i:i+p)'+sigma)*x(i:i+p)*g*e(i+p);
% end
% E=E+e.^2;
% end
%  E=E/loop;
%  subplot(5,1,5);
%  plot(y);  grid;
%  xlabel('����n');
%  ylabel('��ֵ');
%  title('����Ӧ�˲�����ź�');
% 
%  figure;
%  plot(E(p+1:N));  grid;
%  
%  figure;
%  plot(d(p+1:N),'g');  grid;   
%  hold on;
%  plot(y(p+1:N),'b');   grid;
%  hold on;
%  ds=y-d;
%  plot(ds(p+1:N),'r');  grid;
%  legend('�����ź�','����Ӧ�˲�������ź�','���ֵ');
% 
% 
% %--------IPNLMS-------------%
% close all;
% clear  all;
% clc;
% p=9;
% loop=100;
% N=500;
% n=0:N-1;
% x=sin(600.1*pi*n);
% x1=0.8*sin(600.1*pi*n);
% s=1/6*sin(2000.1*pi*n);
% d=s+x1;
% 
% figure;
% subplot(5,1,1);
% plot(x);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('Զ�������ź�');
% 
% subplot(5,1,2);
% plot(x1);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% subplot(5,1,3);
% plot(s);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('���������ź�');
% 
% subplot(5,1,4);
% plot(d);  grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% y=zeros(size(n));
% w=zeros(1,p+1);
% e=zeros(size(n));
% r=zeros(1,p+1);
% E=zeros(size(n));
% u=1.0;          sigma=1.0e-12;  
% ro=1/(p+1);  sigmap=0.01;
% ipsino=-1;  %(����NLMS)-1<ipsino<1(����IPNLMS)
% for q=1:loop
% for i=1:N-p
%     y(i+p)=w*x(i:i+p)';
%     e(i+p)=d(i+p)-y(i+p);
%     temp=abs(w);
%     temp=max(temp);
%     rmax=ro*max(temp,sigmap);
%     for k=1:p+1
%         r(k)=max(rmax,w(k));
%     end
%     g=(p+1)*r./sum(abs(r));
%     for k=1:p+1
%         g(k)=(1-ipsino)/2+(1+ipsino)/2*g(k);
%     end
%     g=diag(g);
%     w=w+u/(x(i:i+p)*g*x(i:i+p)'+sigma)*x(i:i+p)*g*e(i+p);
% end
% E=E+e.^2;
% end
%  E=E/loop;
%  subplot(5,1,5);
%  plot(y);  grid;
%  xlabel('����n');
%  ylabel('��ֵ');
%  title('����Ӧ�˲�����ź�');
% 
%  
%  figure;
%  plot(E(p+1:N));  grid;
%  
%  figure;
%  plot(d(p+1:N),'g');  grid;   
%  hold on;
%  plot(y(p+1:N),'b');   grid;
%  hold on;
%  ds=y-d;
%  plot(ds(p+1:N),'r');  grid;
%  legend('�����ź�','����Ӧ�˲�������ź�','���ֵ');
% 
% 
% 
% %--------------------UPNLMS-------------------%
% close all;
% clear  all;
% clc;
% p=9;
% loop=100;
% N=300;
% n=0:N-1;
% x=sin(600.1*pi*n);
% x1=0.8*sin(600.1*pi*n);
% s=1/6*sin(2000.1*pi*n);
% d=s+x1;
% 
% figure;
% subplot(5,1,1);
% plot(x);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('Զ�������ź�');
% 
% subplot(5,1,2);
% plot(x1);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% subplot(5,1,3);
% plot(s);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('���������ź�');
% 
% subplot(5,1,4);
% plot(d);  grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% y=zeros(size(n));
% w=zeros(1,p+1);
% e=zeros(size(n));
% r=zeros(1,p+1);
% P=zeros(1,p+1);
% E1=zeros(size(n));
% E2=zeros(size(n));
% ERLE=zeros(size(n));
% u=1.0;       sigma=1.0e-12;  beita=1.0e-12;
% ro=1/(p+1);  sigmap=0.01;
% ipsino=1;  %(����NLMS)-1<ipsino<1(����IPNLMS)
% for q=1:loop
%   for i=1:N-p
%     y(i+p)=w*x(i:i+p)';
%     e(i+p)=d(i+p)-y(i+p);
%     temp=abs(w);
%     temp=max(temp);
%     rmax=ro*max(temp,sigmap);
%     for k=1:p+1
%         r(k)=max(rmax,w(k));
%     end
%     g=(p+1)*r./sum(abs(r));
%     for k=1:p+1
%         g(k)=(1-ipsino)/2+(1+ipsino)/2*g(k);
%     end
%     g=diag(g);
%     P=w+u/(x(i:i+p)*g*x(i:i+p)'+sigma)*e(i+p)*x(i:i+p)*g*g; 
%     temp=d(i+p)-P*x(i:i+p)';
%     w=P+u/(beita+x(i:i+p)*x(i:i+p)')*x(i:i+p)*temp;
%   end
%    E1=E1+e.^2;  E2=E2+e;
% end
% 
% E1=E1/loop;   %������ֵ
% E2=E2/loop   %����ֵ
%  up=0;   
%  down=0;
% for k=p+2:N
%     for i=1:p+1
%         up=up+d(k-i)^2;
%         down=down+E2(k-i)^2;
%     end
%   ERLE(k)=10*log10(up/down);
%     %   ERLE(k)=up/down;
% end
%  plot(y);  grid;
%  xlabel('����n');
%  ylabel('��ֵ');
%  title('����Ӧ�˲�����ź�');
% 
%  figure;
%  plot(E1(p+1:N));  grid;
%  figure;
% plot(ERLE(p+2:N));
%  figure;
%  plot(d(p+1:N),'g');  grid;   
%  hold on;
%  plot(y(p+1:N),'b');   grid;
%  hold on;
%  ds=y-d;
%  plot(ds(p+1:N),'r');  grid;
%  legend('�����ź�','����Ӧ�˲�������ź�','���ֵ');
% 
% 
% %----------------------RLS----------------------------%
% %-----------RLS�㷨------------%
% close all;
% clear  all;
% clc;
% p=9;
% M=p+1;
% loop=5;
% N=500;
% n=0:N-1;
% x=sin(600.1*pi*n);
% x1=0.8*sin(600.1*pi*n);
% s=1/6*sin(2000.1*pi*n);
% d=s+x1;
% 
% figure;
% subplot(5,1,1);
% plot(x);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('Զ�������ź�');
% 
% subplot(5,1,2);
% plot(x1);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% subplot(5,1,3);
% plot(s);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('���������ź�');
% 
% subplot(5,1,4);
% plot(d);  grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% Lambda = 0.998 ;    % ��������
% Delta = 0.001 ; 
% y=zeros(size(n));
% w=zeros(1,p+1);
% e=zeros(size(n));
% E=zeros(size(n));
% P = Delta * eye ( M, M ) ;
% 
% for q=1:loop
%    for i=1:N-p
%     % y(i+p)=w*x(i:i+p)';
%     % e(i+p)=d(i+p)-y(i+p);
%     pi_=x(i:i+p)*P;
%     k=Lambda + pi_ *x(i:i+p)';
%     K = pi_'/k; 
%     y(i+p)=w*x(i:i+p)';
%     e(i+p)=d(i+p)-y(i+p);
%     w = w + K' * e(i+p);
%     PPrime = K * pi_ ;
%     P = ( P - PPrime ) / Lambda ;
%    end
%    E=E+e.^2;
% end
%  E=E/loop;
%  subplot(5,1,5);
%  plot(y);  grid;
%  xlabel('����n');
%  ylabel('��ֵ');
%  title('����Ӧ�˲�����ź�');
% 
%  figure;
%  plot(E);  grid;
%  figure;
%  plot(d(p+1:N),'g');  grid;   
%  hold on;
%  plot(y(p+1:N),'b');   grid;
%  hold on;
%  ds=y-d;
%  plot(ds(p+1:N),'r');  grid;
%  legend('�����ź�','����Ӧ�˲�������ź�','���ֵ');
% 
% 
% 
% %------------------------------APA-------------------------------%
% %-----------APA�㷨������------------%
% close all;
% clear  all;
% clc;
% p=9;
% K=10;
% loop=5;
% N=100;
% n=0:N-1;
% x=sin(600.1*pi*n);
% x1=0.8*sin(600.1*pi*n);
% s=1/6*sin(2000.1*pi*n);
% d=s+x1;
% 
% figure;
% subplot(5,1,1);
% plot(x);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('Զ�������ź�');
% 
% subplot(5,1,2);
% plot(x1);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% subplot(5,1,3);
% plot(s);   grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('���������ź�');
% 
% subplot(5,1,4);
% plot(d);  grid;
% xlabel('����n');
% ylabel('��ֵ');
% title('�����ź�');
% 
% w0=zeros(p+1,1);
% mu=0.1;
% epsilon=1.0e-12;
% [y,e,w]=apa_func(x',d',p+1,K,w0,mu,epsilon);
% e=e.^2;
%  subplot(5,1,5);
%  plot(y);  grid;
%  xlabel('����n');
%  ylabel('��ֵ');
%  title('����Ӧ�˲�����ź�');
% 
%  figure;
%  plot(e);  grid;
%  figure;
%  plot(d(p+1:N),'g');  grid;   
%  hold on;
%  plot(y(p+1:N),'b');   grid;
%  hold on;
% ds=y-d';
% plot(ds(p+1:N),'r');  grid;
%  legend('�����ź�','����Ӧ�˲�������ź�','���ֵ');
% 

%-----------------FUNCTION����------------------------------------%

function [y,e,w] = apa_func(x,d,N,K,w0,mu,epsilon)
% ----------------
% input parameters
% ----------------
% x : Lx1 input signal
% d : Lx1 desired response 
% N : filter length
% K : APA order
% w0 : Nx1 initialization
% mu : step-size parameter
% epsilon: regularization parameter (to avoid divide by zero)
% ----------------
% function outputs
% ----------------
% e : Lx1 coefficient output error vector
% w : LxN adaptive filter coefficients

L = length(x);
w = zeros(L,N);
e = zeros(L,1);
y = zeros(L,1);

w(1,:) = w0';
xvec = zeros(N,1);
X = zeros(K,N);
dvec = zeros(K,1);

for i = 1:L-1
xvec = [x(i);xvec(1:N-1)];
X = [xvec';X(1:K-1,:)];
y(i) =w(i,:)*xvec;
e(i) = d(i)-y(i);
dvec = [d(i);dvec(1:K-1)];
evec = dvec - X*w(i,:)';
upd = mu*X'*inv(epsilon*eye(K)+X*X')*evec;
w(i+1,:) = w(i,:) + upd';
end
xvec = [x(L);xvec(1:N-1)];
y(L)=w(L,:)*xvec;
e(L) = d(L)-y(L);
end