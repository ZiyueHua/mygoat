function [dMdt] = GOATEvolution(t,M,H0,Hs,A,W,P,n_har)
% This function computes unitary evolution from shrodinger
% equation
% persistent itr;
% if isempty(itr)
%    itr=0;
% end
% itr = itr + 1;
% 
% if mod(itr,1000) == 0
%     fprintf('perevol time: %f\n', toc);
%     fprintf('GOATEvolution Iteration %d\n', itr);
% end

f=zeros(length(Hs),1);
f1=zeros(length(Hs),1);
f2=zeros(length(Hs),1);

dMdt = zeros((1+length(Hs)*n_har+length(Hs))*int32(length(H0)^2),1);

for k = 1 : length(Hs)
    for m = 1 : n_har
        f(k) = f(k) + A(k,m) * sin( m * t.* W(k) + P(k) );
        f1(k) = f1(k)  + A(k,m) * m * t * cos( m * t * W(k) + P(k));
        f2(k) = f2(k)  + A(k,m) * cos( m * t * W(k) + P(k));
    end     
end

H = H0;
for i = 1 : length(Hs)
    H = H + f(i) * Hs{i};
end
     
% Lt = kron(eye(4),H);

U = reshape( M(1:int32(length(H0)^2),1), size(H0) );
U = -1i*H*U;
dMdt(1:int32(length(H0)^2),1) = U(:);


%% This part related to derivatives with respect to waight matrices
for k = 1 : length(Hs)
    for m = 1 : n_har
        pos = ( ( k - 1 ) * n_har + m ) * int32(length(H0)^2);
        dUdtpkm = reshape( M( pos + 1 : pos + int32(length(H0)^2), 1 ), size(H0) );
        dUdtpkm = -1i * ( sin( m * t * W(k) + P(k) ) * Hs{k} * U + H * dUdtpkm );
        dMdt( pos + 1 : pos + int32(length(H0)^2) , 1 ) = dUdtpkm(:);    
    end
end

%% This part related to derivatives with respect to fourier base frequencies
%  for each control
for k = 1 : length(Hs)
    pos = ( length(Hs) * n_har + k ) * int32(length(H0)^2);
    dUdtpk = reshape( M( pos + 1 : pos + int32(length(H0)^2), 1 ), size(H0) );
    dUdtpk =  -1i * ( f1(k) * Hs{k} * U + H * dUdtpk );
    dMdt( pos + 1 : pos + int32(length(H0)^2), 1 ) = dUdtpk(:);
end

for k = 1 : length(Hs)
    pos = ( length(Hs) * n_har + length(Hs) + k ) * int32(length(H0)^2);
    dUdtpk = reshape( M( pos + 1 : pos + int32(length(H0)^2), 1 ), size(H0) );
    dUdtpk = -1i * ( f2(k) * Hs{k} * U + H * dUdtpk );
    dMdt( pos + 1 : pos + int32(length(H0)^2), 1 ) = dUdtpk(:);
end

end


