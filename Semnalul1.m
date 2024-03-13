t=Rogoz(:,1)
u=Rogoz(:,2)
y1=Rogoz(:,3)

plot(t,[u,y1])%reprezentarea grafica a primului set de date

%% Identificarea cu ajutorul fenomenului de rezonanta
 i1=316; %ymax
 i2=329; %ymin
 i3=309; %umax
 i4=323; %umin
 k=mean(y1)/mean(u)
 Mr=(y1(i1)-y1(i2))/(u(i3)-u(i4))/k
 %pentru semnalul de iesire
 Tr1=2*(t(i2)-t(i1))
 %pentru semnalul de intrare

Tr2=2*(t(i4)-t(i3))
wr2=2*pi/Tr2
zetta2=sqrt((Mr-sqrt(Mr^2-1))/2/Mr)
wn2=wr2/sqrt(1-2*zetta2^2)

num2=k*wn2^2
den2=[1 2*zetta2*wn2 wn2^2]
A=[0,1;-wn2^2,-2*zetta2*wn2]
B=[0;k*wn2^2]
C=[1,0]
D=[0]
sys=ss(A,B,C,D);
ysim=lsim(sys,u,t,[y1(1),(y1(2)-y1(1))/(t(2)-t(1))]);

figure 
plot(t,y1,t,ysim)
J=norm(y1-ysim)/sqrt(length(y1)) %eroare medie patratica
eMPN=norm(y1-ysim)/norm(y1-mean(y1))*100%eroare medie patratica normalizata

% figure 
% nyquist(num2,den2)
%% Identificare Nyquist
 i5=343;
 i6=356;
 i7=337;
 i8=349;

dt=t(i7)-t(i5)
Tn=2*(t(i8)-t(i7))
wn3=2*pi/Tn
ph=dt*wn3*180/pi 


M1=(y1(i5)-y1(i6))/(u(i7)-u(i8))/k 
zetta3=k/M1/2
num3=k*wn3^2
den3=[1 2*zetta3*wn3 wn3^2 ]
A=[0,1;-wn3^2,-2*zetta3*wn3] 
B=[0;k*wn3^2]
C=[1,0]
D=[0]
sys=ss(A,B,C,D);
ysim=lsim(sys,u,t,[y1(1),(y1(2)-y1(1))/(t(i2)-t(i1))]);
figure 
plot(t,y1,t,ysim)
J=norm(y1-ysim)/sqrt(length(y1)) %eroare medie patratica
eMPN=norm(y1-ysim)/norm(y1-mean(y1))*100%eroarea medie patratica normalizata

figure
nyquist(num3,den3)
%%
dt = t(2)-t(1) %timp de achizitie, diferenta dintre 2 momente de timp consecutive
d_id = iddata(y1,u,dt) %datele de identificare;semnal iesire, intrare, pas de achizitie
%% Model calculat cu arx
Marx = arx(d_id,[2,2,1])%nA,nB,nk
figure;resid(d_id, Marx, 4) 
figure;compare(d_id, Marx) 

%% a doua metoda de validare - cu variab instrumentale MVI
Mvi = iv4(d_id, [2,2,1])%nA,nB,nk
figure;resid(d_id, Mvi,4)%validare statica
figure;compare(d_id,Mvi)%grad de suprapunere
Hz = tf(Mvi.B, Mvi.A, dt) %functie de transfer in z sau modelul in discret;
Hs = d2c(Hz, 'zoh')%functia de transfer in continu

%% Model cu output error
Moe = oe(d_id,[2,2,1]) %nb,nf,nk 
figure;resid(d_id,Moe,5)%validare statica
figure;compare(d_id,M_oe) %grad de suprapunere

%% model cu armax
Marmax = armax(d_id,[2,2,2,1]) %na, nb,nc, nd(dimensiunea ferestrei alunecatoare)
figure;resid(d_id,Marmax,5)%validare statica
figure;compare(d_id,Marmax)%grad de suprapunere
Hz = tf(Marmax.B, Marmax.A, dt) %functie de transfer in z sau modelul in discret;
Hs = d2c(Hz, 'zoh')%functie de transfer in continu

%% 
%Mpem = pem(d_id, 2) 
Marx_pem=pem(Marx,d_id)
figure
resid(d_id, Marx_pem)
figure
compare(d_id, Marx,Marmax,Mvi,Moe,Marx_pem)
%% metoda de minimizare a err de pred obtinut cu n4sid
n4sid(d_id, 1:10)
Mn4sid=n4sid(d_id,2)
figure
resid(d_id, Mn4sid)
figure
compare(d_id, Marx,Marmax,Mvi,Moe,Marx_pem,Mn4sid)
%% Modelul cu arx a avut cea mai mare eroare

m_arx_pem= pem(Marx,d_id)%rafinarea prin utilizarea functiei pem
figure; resid(d_id,m_arx_pem)%validarea statistica
figure;compare(d_id,m_arx_pem)%gradul de suprapunere
