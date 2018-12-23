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
errs = sum(xor(symbols, symbols2));
shiftedErrs = sum(xor(symbols, symbols3));
    
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


