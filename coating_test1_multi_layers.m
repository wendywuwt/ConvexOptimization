clear
clc
layer=10;
n=layer/2+1;
omega=linspace(0.5*pi,pi,1024);
na=2;
np=1;
C=(na-np)/(na+np);
for j=1:1024
    A(j,:)=g2(omega(j),n);
end
cvx_begin
    variable rho(n)
    minimize(max(abs(A*rho)))
%     rho1=rho;
%     for j=layer/2+2:layer+1
%         rho1(j)=rho(layer+2-j);
%     end

    subject to
        g2(0,n)*rho==C

        abs(rho)<=ones(n,1)

cvx_end

omega1=linspace(0,pi,2048);
for j=1:2048
    B(j,:)=g2(omega1(j),n);
end
gama=(abs(B*rho));
for j=layer/2+2:layer+1
    rho(j)=rho(layer+2-j);
end
refr(1)=na;
for j=2:layer+1
    refr(j)=(1-rho(j-1))/(1+rho(j-1))*refr(j-1);
end
refr(layer+2)=np;
for j=2049:4096
    gama(j)=gama(4096+1-j);
end
pointer=1;
for j=2:length(refr)-1
       if rho>0
        for k=pointer:pointer+100/refr(j)
            thickness(k)=refr(j);
        end
    else
        for k=pointer:pointer+100/refr(j-1)
            thickness(k)=refr(j-1);
        end
        pointer=k+1;
        for k=pointer:pointer+100/refr(j)
            thickness(k)=refr(j);
        end
    end
    pointer=k+1;
end
% for k=pointer:pointer+100/refr(j)
%             thickness(k)=refr(j);
%         end
x_for_thick=linspace(0,length(thickness)/100/4,length(thickness));
% U=[1 0;0 1];
% for j=1:2048
%     z(j)=exp(2*i*omega1(j));
%     for k=1:length(rho)
%         U=z(j)/(rho(k)+1)*[1 rho(k)*z(j)^-1;rho(k) z(j)^-1]*U;
%     end
%     gama1(j)=abs(1/U(1,1));
% end
figure(1)
[ax,h1,h2]=plotyy(linspace(0,2,4096),20*log10(gama),linspace(0,2,4096),gama*100);
set(ax(1),'YLim',[-150,0])
set(ax(2),'YLim',[-10,100])
xlabel('\itf / f_{\rm0}','fontsize',15,'fontname','times new roman')
ylabel(ax(1),'|\Gamma( \itf \rm )|^2/dB','fontsize',15,'fontname','times new roman')
ylabel(ax(2),'|\Gamma( \itf \rm )|^2/%','fontsize',15,'fontname','times new roman')
title([num2str(layer),' Layers Reflection Coefficient',' n_a=',num2str(na),' n_b=',num2str(np)],'fontsize',15,'fontname','times new roman')
figure(2)
area(x_for_thick,thickness)
xlabel('Thickness/( \it\lambda / \lambda_{\rm0} \rm)','fontsize',15,'fontname','times new roman')
ylabel('Refraction Factor \itn','fontsize',15,'fontname','times new roman')
title([num2str(layer),' Strcuture of Coating',' n_a=',num2str(na),' n_b=',num2str(np)],'fontsize',15,'fontname','times new roman')

