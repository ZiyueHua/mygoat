function [Fidelity,g,res]=GOATPulse(u,H0,Hs,inis,trgs,tspan,n_har,h,penaltyl,filename)
save( filename, 'u')

persistent itr;
if isempty(itr)
   itr=0;
end
itr = itr + 1;
tic
fprintf('GOATPulse1 Iteration %d\n', itr);

phi = zeros(length(trgs),1);
g = zeros(n_har*length(Hs) + 2 * length(Hs), 1);

for i = 1 : length(inis)
    U=eye(length(H0));
    M0=zeros((1+length(Hs)*n_har+2*length(Hs))*int32(length(U)^2),1);
    M0(1:int32(length(U)^2),1)=U(:);
    A=reshape(u(1:length(Hs)*n_har,1), length(Hs), n_har);
    W=u(length(Hs)*n_har+1:length(Hs)*n_har+length(Hs),1);
    P=u(length(Hs)*n_har+length(Hs)+1:end,1);

    %% solve ode
    clear GOATEvolution;
    opt = odeset('RelTol',1e-12,'AbsTol',1e-12,'Stats','off');
    [~,M] = ode45(@(t,M) GOATEvolution(t,M,H0,Hs,A,W,P,n_har), tspan, M0, opt);

    %% calculate result

    U = reshape( M( length(tspan), 1:int32(length(U)^2) ), size(H0) ); 

    phi(i,1) = abs(trgs{i}'* U * inis{i})^2;

    for k = 1 : length(Hs)
        for m = 1 : n_har
            pos = ( ( k - 1 ) * n_har + m ) * int32(length(U)^2);
            dUpkm = reshape( M( length(tspan), pos + 1 : pos + int32(length(U)^2) ), size(H0) );
            phi0pkm = - 2 * real( trgs{i}' * dUpkm * inis{i});
            g(( k - 1 ) * n_har + m) = g(( k - 1 ) * n_har + m) + phi0pkm;
        end
    end

    for k = 1 : length(Hs)
        pos = ( length(Hs) * n_har + k ) * int32(length(U)^2);
        dUpk = reshape( M( length(tspan), pos + 1 : pos + int32(length(U)^2) ), size(H0) );
        phi0pk = - 2 * real(trgs{i}' * dUpk * inis{i});
        g(length(Hs) * n_har + k) = g(length(Hs) * n_har + k) + phi0pk;
    end
    
    for k = 1 : length(Hs)
        pos = ( length(Hs) * n_har + length(Hs) + k ) * int32(length(U)^2);
        dUpk = reshape( M( length(tspan), pos + 1 : pos + int32(length(U)^2) ), size(H0) );
        phi0pk = - 2 * real(trgs{i}' * dUpk * inis{i});
        g(length(Hs) * n_har + length(Hs) + k) = g(length(Hs) * n_har + length(Hs) + k) + phi0pk;
    end

end


phi0 = 0;
for i = 1 : length(inis)
    phi0 = phi0 + abs(phi(i,1))/length(inis);
end

Phkm=exp((u(1:n_har*length(Hs), 1)./h(1:n_har*length(Hs), 1)).^2)-1;
dPhkm=2*u(1:n_har*length(Hs), 1)./h(1:n_har*length(Hs), 1) .* (Phkm+1);
Phk=exp((u(n_har*length(Hs)+1:n_har*length(Hs)+length(Hs), 1)./h(n_har*length(Hs)+1:end, 1)).^2)-1;
dPhk=2*u(n_har*length(Hs)+1:n_har*length(Hs)+length(Hs), 1)./h(n_har*length(Hs)+1:end, 1) .* (Phk+1);
dP = cat(1, dPhkm, dPhk, zeros(length(Hs),1));

Fidelity=1-phi0+penaltyl*(sum(sum(Phkm))+sum(Phk));

if nargout>1
    g = g./length(inis) - penaltyl .* dP;
end
if nargout>2
    res=U*inis{1};
end

fprintf('total time: %f\n', toc);

end