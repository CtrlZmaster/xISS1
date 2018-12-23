% Projekt: ISS Projekt 2018/2019
% Súbor: moje_reseni.m
% Autor: Michal Pospíšil (xpospi95@stud.fit.vutbr.cz)

% Loading signal package in Octave
pkg load signal

%% Úloha 1
disp("Úloha 1");
% Výpis informácii o súbore 
info = audioinfo("xpospi95.wav");
N = info.TotalSamples;
Fs = info.SampleRate;
disp("Počet symbolov:");
disp(N/16); % Zadanie: 1 symbol reprezentuje 16 vzoriek

disp("Vzorkovacia frekvencia:");
disp(Fs);

disp("Trvanie:");
disp(info.Duration);

% Načítanie signálu s
s = audioread("xpospi95.wav");

%% Úloha 2
disp("Úloha 2");
% Priemerovanie hodnôt z každých 16 vzoriek
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
figure('Name','Úloha 2');
plot(linspace(0,0.020,320), s(1:320), "b");
title("Signál a dekódované symboly"); xlabel("t"); grid;

hold on 
stem(linspace(0.0005,0.0195,20), symbols(1:20), "g");
leg1 = legend("s[n]","Symboly");
legend(leg1, 'Location','southoutside','Orientation','horizontal');
hold off

% Overenie zhody so súborom
fileID = fopen('xpospi95.txt','r');
formatSpec = '%d';
file_symbols = fscanf(fileID,formatSpec);

if(sum(xor(symbols, file_symbols)) > 0)
  disp('Symboly boli dekódované správne.');
else
  disp('Symboly NEBOLI dekódované správne.');
endif

%% Úloha 3
disp("Úloha 3");
% Koeficienty filtra
B = [0.0192 -0.0185 -0.0185 0.0192];
A = [1.0000 -2.8870 2.7997 -0.9113];
figure('Name','Úloha 3');
zplane(B,A); title("Nuly a póly filtra"); xlabel('Re'); ylabel('Im');

% Stabilita filtra
p = roots(A); 
if (isempty(p) | abs(p) < 1) 
  disp('Filter je stabilný.')
else
  disp('Filter je NESTABILNÝ.')
endif

%% Úloha 4
disp("Úloha 4");
% Frekvenčná charakteristika
[H, w] = freqz(B,A,Fs); f = (0:Fs-1)/ Fs * Fs / 2;
figure('Name','Úloha 3 - Frekvenčná charakteristika filtra');
subplot(2,1,1); plot(f,abs(H)); title("Modul kmitočtovej charakteristiky filtra"); xlabel('f'); ylabel('|H(f)|'); grid;
subplot(2,1,2); plot(f,angle(H)); title("Faza kmitočtovej charakteristiky filtra"); xlabel('f'); ylabel('arg H(f)'); grid;

%% Úloha 5
disp("Úloha 5");
% Filtrovanie signálu
ss = filter(B,A,s);

figure('Name','Úloha 5');
plot(linspace(0,0.020,320), s(1:320), "b", linspace(0,0.020,320), ss(1:320), "g");
title("Rozdiel signálov po filtrovaní"); xlabel('t'); grid;
leg2 = legend("s[n]","ss[n]");
legend(leg2, 'Location','southoutside','Orientation','horizontal');

% Výpočet korelačného koeficientu 
[cor,lag] = xcorr(ss,s);

% Výber najvyššieho koeficientu, tam sa signály najviac zhodujú - sú zosynchronizované - z MATLAB pomocnej stránky
[~,I] = max(abs(cor));
lagDiff = lag(I); 
timeDiff = lagDiff/Fs;
disp("Oneskorenie signálu ss [s]:")
disp(timeDiff);

%% Úloha 6
disp("Úloha 6");
% Zarovnanie filtrovaného signálu
ssa = ss(lagDiff+1:N-lagDiff);

% Vykreslenie všetkých signálov
figure('Name','Úloha 6');
plot(linspace(0,0.020,320), s(1:320), "b", linspace(0,0.020,320), ss(1:320), "g", linspace(0,0.020,320), ssa(1:320), "m");
title("Posunutie, filtrovanie a zarovnanie"); xlabel("t"); grid;
legend('Location','southoutside','Orientation','horizontal');
% Priemerovanie hodnôt z každých 16 vzoriek
i = 1;
j = 1;
shiftedthresh = mean(ssa);
disp("Modifikovaný prah:");
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
stem(linspace(0.0005,0.0195,20), symbols2(1:20), "r");
leg3 = legend("s[n]","ss_{shifted}[n]","Zarovnaný ss_{shifted}[n]","Symboly zo zarovnaného ss_{shifted}[n]");
legend(leg3, 'Location','southoutside','Orientation','horizontal');
hold off

%% Úloha 7
disp("Úloha 7");
% Rátanie počtu nezhôd
errs = sum(xor(symbols(1:end-1), symbols2));
shiftedErrs = sum(xor(symbols(1:end-1), symbols3));
    
disp("Chyba pri dekódovaní z filtrovaného signálu [%]:");
disp(errs/(N/16)*100);

disp("Chyba pri dekódovaní z filtrovaného signálu (posunutý prah) [%]:");
disp(shiftedErrs/(N/16)*100);

%% Úloha 8
disp("Úloha 8");

% Graf spektier zadaného a posunutého signálu
figure('Name','Úloha 8');
fftRes1 = abs(fft(s));
fftRes1 = fftRes1(1:Fs/2);
subplot(2,1,1); plot(0:Fs/2-1, fftRes1);
title("Spektrum signálu s[n]"); xlabel("f"); ylabel("s[n]"); grid;

fftRes2 = abs(fft(ss));
fftRes2 = fftRes2(1:Fs/2);
subplot(2,1,2); plot(0:Fs/2-1, fftRes2);
title("Spektrum signálu ss[n]"); xlabel("f"); ylabel("s[n]"); grid;


%% Úloha 9
disp("Úloha 9");
figure('Name','Úloha 9');

hist(s,15,1);
title("Odhad PDF s[n]"); xlabel("x"); ylabel("p(x z intervalu, n)"); grid;
%xticks(-1:0.1:1);
set(gca,'XTick',-1:0.1:1); % LINUX nepozna xticks

%% Úloha 10
disp('Úloha 10');

[cor10, lag10] = xcorr(s, 50, 'biased');

figure('Name','Úloha 10');
plot(lag10, cor10);
title("Korelačné koeficienty"); xlabel("k"); ylabel("R[k]"); grid;

%% Úloha 11
disp('Úloha 11');

disp('R[0]:');
disp(cor10(51));

disp('R[1]:');
disp(cor10(52));

disp('R[16]:');
disp(cor10(67));

%% Úloha 12
disp('Úloha 12');

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

figure('Name','Úloha 12');
imagesc(chlieviky, chlieviky, p12);
title("Časový odhad združenej funkcie hustoty rozdelenia pravdepodobnosti");
colorbar;

%% Úloha 13
disp('Uloha 13');
disp("Dôkaz urobila funkcia hist2 a je vytlačený v predchádzajucej úlohe.");

%% Úloha 14
disp('Úloha 14');
disp('R[1]:');
disp(r12);


%% Uloženie obrázkov do súboru
for i=1:9
clear name;
name = ["fig" num2str(i) ".svg"];
saveas (i, name);
endfor

for i=1:9
clear name;
name = ["fig" num2str(i) ".png"];
saveas (i, name);
endfor