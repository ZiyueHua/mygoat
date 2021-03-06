function [dMdt] = GOATEvolution(t,M,H0,Hs,A,W,P,n_har)
% This function computes unitary evolution from shrodinger
% equation
% persistent itr;
% if isempty(itr)
%    itr=0;
% end
% itr = itr + 1;
% 
% if mod(itr,100) == 0
%     fprintf('perevol time: %f\n', toc);
%     fprintf('GOATEvolution Iteration %d\n', itr);
% end

f=zeros(length(Hs),1);
f1=zeros(length(Hs),1);
f2=zeros(length(Hs),1);

dMdt = zeros((1+length(Hs)*n_har+length(Hs))*length(H0),1);

for k = 1 : length(Hs)
    for m = 1 : n_har
        f(k) = f(k) + A(k,m) * sin( m * t.* W(k) + P(k) );
    end     
end

H = H0;
for i = 1 : length(Hs)
    H = H + f(i) * Hs{i};
end
     
% Lt = kron(eye(4),H);

Uini = M(1:length(H0),1);
dMdt(1:length(H0),1) = -1j*H*Uini;


%% This part related to derivatives with respect to waight matrices
for k = 1 : length(Hs)
    for m = 1 : n_har
        f1(k) = f1(k)  + A(k,m) * m * t * cos( m * t * W(k) + P(k));
        f2(k) = f2(k)  + A(k,m) * cos( m * t * W(k) + P(k));
        pos = ( ( k - 1 ) * n_har + m ) * length(H0);
        dUdtpkmini = M( pos + 1 : pos + length(H0), 1 );
        dMdt( pos + 1 : pos + length(H0) , 1 ) = -1i * ( sin( m * t * W(k) + P(k) ) * Hs{k} * Uini + H * dUdtpkmini );    
    end
end

%% This part related to derivatives with respect to fourier base frequencies
%  for each control
for k = 1 : length(Hs)
    pos = ( length(Hs) * n_har + k ) * length(H0);
    dUdtpkini = M( pos + 1 : pos + length(H0), 1 );    
    dMdt( pos + 1 : pos + length(H0), 1 ) = -1i * ( f1(k) * Hs{k} * Uini + H * dUdtpkini );
end

for k = 1 : length(Hs)
    pos = ( length(Hs) * n_har + length(Hs) + k ) * length(H0);
    dUdtpkini = M( pos + 1 : pos + length(H0), 1 );    
    dMdt( pos + 1 : pos + length(H0), 1 ) = -1i * ( f2(k) * Hs{k} * Uini + H * dUdtpkini );
end

end


