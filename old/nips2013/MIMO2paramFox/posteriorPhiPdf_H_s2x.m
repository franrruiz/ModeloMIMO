function [pdf grad] = posteriorPhiPdf_H_s2x(H,X,ZR,chi,tau,nu)
% Z: Tx(MQ)
% Phid: (MQ)x1
% Xd = Tx1

% Devuelve "menos la log-pdf" (no normalizada) del posterior de las Phi's,
% ademas del gradiente
D=size(X,2);
[T MK] = size(ZR);

pdf = (tau+MK*D/2)*log(nu+chi/2*trace(H'*H))+(tau+T*D/2)*log(nu+1/2*trace((X-ZR*H)'*(X-ZR*H)));

grad = (tau+MK*D/2)*chi*H/(nu+chi/2*trace(H'*H))+(tau+T*D/2)*(ZR'*ZR*H-ZR'*X)/(nu+1/2*trace((X-ZR*H)'*(X-ZR*H)));
