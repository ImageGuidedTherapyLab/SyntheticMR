
function[T2]=qalasT2calc(M0,T1,Mmeas,TR,TE_T2prep,flipAngle,nacq,dt)

Mend=M0-(M0-Mmeas(end)./sind(flipAngle)).*exp(-dt(end)./T1);
T2=-TE_T2prep./log((Mmeas(1)./sind(flipAngle))./Mend);

% 
% star=(1-exp(-TR./T1))./(1-cosd(flipAngle)*exp(-TR./T1));
% 
% % T2 sensitization
% M(3)=Mmeas(1);
% M(4)=M0-(M0-M(3)).*exp(-dt(4)./T1);
% 
% % T1 sensitization
% M(5)=-M(4);             % Assume perfect inversion? 100 ms inversion time
% M(6)=M0-(M0-M(5))*exp(-dt(6)/T1);
% 
% % Post T1 sens acquisitions
% for iii=1:nacq-1
%     M(5+2*iii)=M0.*star-(M0.*star-M(4+2*iii)).*exp(-dt(5+2*iii)./(T1.*star));
%     M(6+2*iii)=M0-(M0-M(5+2*iii)).*exp(-dt(6+2*iii)./T1);
% end
% 
% M(2)=M0.*star-(M0.*star-Mmeas(1))./exp(-dt(3)./(T1.*star));
% T2=-TE_T2prep./log(M(2)./M(end));

end