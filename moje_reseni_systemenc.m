% Projekt: ISS Projekt 2018/2019
% S�bor: moje_reseni.m
% Autor: Michal Posp�il (xpospi95@stud.fit.vutbr.cz)

% Loading signal package in Octave
pkg load signal

%% �loha 1
disp("Uloha 1");
% V�pis inform�cii o s�bore 
info = audioinfo("xpospi95.wav");
N = info.TotalSamples;
Fs = info.SampleRate;
disp("Pocet symbolov:");
disp(N/16); % Zadanie: 1 symbol reprezentuje 16 vzoriek

disp("Vzorkovacia frekvencia:");
disp(Fs);

disp("Trvanie:");
disp(info.Duration);

% Na��tanie sign�lu s
s = audioread("xpospi95.wav");

%% �loha 2
disp("Uloha 2");
% Priemerovanie hodn�t z ka�d�ch 16 vzoriek
i = 1;
j = 1;
while(i <= N - 16)
  if(mean(s(i:i+16-1)) > 0)
    symbols(j) = 1;
  else
    symbols(j) = 0;
  endif
  i += 16; % Next symbol
  j += 1;
endwhile

% Kreslenie grafu
figure('Name','Uloha 2');
plot(linspace(0,0.020,320), s(1:320), "b");
title("Signal a dekodovane symboly"); xlabel("t"); ylabel("symboly, s[n]"); grid;
hold on 
stem(linspace(0.0005,0.0195,20), symbols(1:20), "g");
hold off

%% �loha 3
disp("Uloha 3");
% Koeficienty filtra
B = [0.0192 -0.0185 -0.0185 0.0192];
A = [1.0000 -2.8870 2.7997 -0.9113];
figure('Name','Uloha 3');
zplane(B,A); title("Nuly a poly filtra"); xlabel('Re'); ylabel('Im');

% Stabilita filtra
p = roots(A); 
if (isempty(p) | abs(p) < 1) 
  disp('Filter je stabilny.')
else
  disp('Filter je nestabilny.')
endif

%% �loha 4
disp("Uloha 4");
% Frekven�n� charakteristika
[H, w] = freqz(B,A,Fs); f = (0:Fs-1)/ Fs * Fs / 2;
figure('Name','Uloha 3 - Frekvencna charakteristika filtra');
subplot(2,1,1); plot(f,abs(H)); title("Modul kmitoctovej char. filtra"); xlabel('f'); ylabel('|H(f)|'); grid;
subplot(2,1,2); plot(f,angle(H)); title("Faza kmitoctovej char. filtra"); xlabel('f'); ylabel('arg H(f)'); grid;

%% �loha 5
disp("Uloha 5");
% Filtrovanie sign�lu
ss = filter(B,A,s);

figure('Name','Uloha 5');
plot(linspace(0,0.020,320), s(1:320), "b", linspace(0,0.020,320), ss(1:320), "g");
title("Rozdiel signalov po filtrovani"); xlabel('t'); grid;

% V�po�et korela�n�ho koeficientu 
[cor,lag] = xcorr(ss,s);

% V�ber najvy��ieho koeficientu, tam sa sign�ly najviac zhoduj� - s� zosynchronizovan� - z MATLAB pomocnej str�nky
[~,I] = max(abs(cor));
lagDiff = lag(I); 
timeDiff = lagDiff/Fs;
disp("Oneskorenie signalu ss [s]:")
disp(timeDiff);

%% �loha 6
disp("Uloha 6");
% Zarovnanie filtrovan�ho sign�lu
ssa = ss(lagDiff+1:N-lagDiff);

% Vykreslenie v�etk�ch sign�lov
figure('Name','Uloha 6');
plot(linspace(0,0.020,320), s(1:320), "b;s[n];", linspace(0,0.020,320), ss(1:320), "g;Filtrovany s[n];", linspace(0,0.020,320), ssa(1:320), "m;Filtrovany a zarovnany s[n];");
title("Posunutie, filtrovanie a zarovnanie"); xlabel("t"); grid;
legend('Location','southoutside','Orientation','horizontal');
% Priemerovanie hodn�t z ka�d�ch 16 vzoriek
i = 1;
j = 1;
shiftedthresh = mean(ssa);
disp("Modifikovany prah:");
disp(shiftedthresh);
while(i <= length(ssa) - 16)
  if(mean(ssa(i:i+16-1)) > 0)
    symbols2(j) = 1;
  else
    symbols2(j) = 0;
  endif
  if(mean(ssa(i:i+16-1)) > shiftedthresh)
    symbols3(j) = 1;
  else
    symbols3(j) = 0;
  endif
  i += 16; % Next symbol
  j += 1;
endwhile

% Doplnenie symbolov do grafu
hold on 
stem(linspace(0.0005,0.0195,20), symbols2(1:20), "r;Symboly z filtrovaneho a posunuteho s[n];");
hold off

%% �loha 7
disp("Uloha 7");
% R�tanie po�tu nezh�d
errs = sum(xor(symbols(1:end-1), symbols2));
shiftedErrs = sum(xor(symbols(1:end-1), symbols3));
    
disp("Chyba pri dekodovani z filtrovanho signalu [%]:");
disp(errs/(N/16)*100);

disp("Chyba pri dekodovani z filtrovanho signalu (posunuty prah) [%]:");
disp(shiftedErrs/(N/16)*100);

%% �loha 8
disp("Uloha 8");

% Graf spektier zadan�ho a posunut�ho sign�lu
figure('Name','Uloha 8');
fftRes1 = abs(fft(s));
fftRes1 = fftRes1(1:Fs/2);
subplot(2,1,1); plot(0:Fs/2-1, fftRes1);
title("Spektrum sign�lu s[n]"); xlabel("f"); ylabel("s[n]"); grid;

fftRes2 = abs(fft(ss));
fftRes2 = fftRes2(1:Fs/2);
subplot(2,1,2); plot(0:Fs/2-1, fftRes2);
title("Spektrum posunut�ho sign�lu s[n]"); xlabel("f"); ylabel("s[n]"); grid;


%% �loha 9
disp("Uloha 9");
figure('Name','Uloha 9');

hist(s,15,1);
title("Odhad PDF s[n]"); xlabel("x"); ylabel("Pocet vzoriek"); grid
xticks(-1:0.1:1);

%% �loha 10
disp('Uloha 10');

[cor10, lag10] = xcorr(s, 50, 'biased');

figure('Name','Uloha 10');
plot(lag10, cor10);
title("Korelacne koeficienty"); xlabel("k"); ylabel("R[k]"); grid;

%% �loha 11
disp('Uloha 11');

disp('R[0]:');
disp(cor10(51));

disp('R[1]:');
disp(cor10(52));

disp('R[16]:');
disp(cor10(67));

%% �loha 12
disp('Uloha 12');

function [h,p,r] = hist2(y1,y2,x); 
% function [h,p,r] = hist2(y1,y2,x); 
%
% y1 and y2 are column vectors with data (must have same dimension)
% x  are centers of histogram bins for the 2 coordinates. 
%    should be equally spaced !!!
% h is the 2-D histogram (y1 distrib in cols, y2 in rows)
% 
% p is the estimate of 2-D prob. dens. function
%
% r is he autocorrelation coefficient computed using the theoretical 
% formula (5-8 in Sebesta) 
%
% ... histogram computation is not very optimized ... 

L = length(x); 
N = length(y1); 

% alloc for hist 
h = zeros(L,L); out = 0; % counter of bad samples...

% make a BIG matrix with all values of x ... will actually do something 
% like vector quantization. 
% first make col vector out of x, then repeat it N times. 
xcol = x(:); bigx = repmat(xcol,1,N); 

% make rows out of y1 , y2 and do the 'quantization' 
yr = y1(:)'; bigy = repmat(yr,L,1);
[dummy,ind1] = min(abs(bigy - bigx)); 
yr = y2(:)'; bigy = repmat(yr,L,1);
[dummy,ind2] = min(abs(bigy - bigx)); 

% now just go through the indices and increment respective cases in h
for ii=1:N,
  d1 = ind1(ii);   d2 = ind2(ii); 
  h(d1,d2) = h(d1,d2) + 1; 
endfor

%%% prob dens: will have to normalize
surf = (x(2) - x(1))^2; % surface of one tile
p = h / N / surf;  

%%%% autocor coeff 
% make col vector out of x and clone it L times. 
x = x(:); X1 = repmat(x,1,L);
% make row vector out of x and clone it L times. 
x = x'; X2 = repmat(x,L,1); 
% now put it together, don't forget to multipl by tile surface
r = sum(sum (X1 .* X2 .* p)) * surf;

%%% check ... 
check = sum(sum (p)) * surf; 
disp(['hist2: check -- 2d integral should be 1 and is ' num2str(check)]); 
endfunction

% Posun o 1 vlavo 
ss12 = s(2:N);

% Orezanie s na rovnaky rozmer
s12 = s(1:N-1);

chlieviky = linspace(-1,1,25);

[h12, p12, r12] = hist2(s12,ss12,chlieviky);

figure('Name','Uloha 12');
imagesc(chlieviky, chlieviky, p12);
colorbar;

%% �loha 13
disp('Uloha 13');
disp("Dokaz urobila funkcia hist2 a je vytlaceny v predchadzajucej ulohe.");

%% �loha 14
disp('Uloha 14');
disp('R[1]:');
disp(r12);