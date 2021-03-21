function [w,corr]=LubichWeights(f,N,h,order)
% [w,corr]=LubichWeights(f,N,h,order)
% f is the Laplace transform of B with B'=b , h is the time step, N is the number of
% weights required, order is the order of the method.

L=12*N;
    function out=finDiff(z) 
        out=0;
        if order<10,
        for jj=1:order
            out=out+(1-z).^jj/jj;
        end
        else
            out=1/2*(z-1./z);
        end
        out=out/h;
    end

    function outt=fLubich(z)
        if order<10,
            outt=f(finDiff(z));
        else
            outt=f(finDiff(z)).*exp(-finDiff(z).*h*2);
        end
    end
            

rho=(1e-7)^(1/L);
j=0:L-1;
fvec=fLubich(rho.*exp(2*1i*pi*j/L));
wvec=fft(fvec)/L;
w=real(wvec(1:N)./rho.^j(1:N));

corr=[];
switch order
    case 2
        corr=-1/2;
    case 3
        corr=[1/12;-7/12];
    case 4
        corr=[-1/24;1/6;-5/8];
end
end



%%% to compute next time step of v(x,t)
% v(:,j+1)=(w(1)*eye(length(xcoord))-A)\(-v(:,1:j)*w(j+1:-1:2)-v(:,end-order+2:end)*corr)

        