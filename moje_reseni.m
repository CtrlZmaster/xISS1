% Projekt: ISS Projekt 2018/2019
% Súbor: moje_reseni.m
% Autor: Michal Pospíšil (xpospi95@stud.fit.vutbr.cz)


% Úloha 1
% Výpis informácii o súbore 
info = audioinfo("xpospi95.wav");
N = info.TotalSamples
disp("Počet symbolov:");
disp(N/16); % Zadanie: 1 symbol reprezentuje 16 vzoriek

disp("Vzorkovacia frekv.:");
disp(info.SampleRate);

disp("Dĺžka:");
disp(info.Duration);

% Načítanie signálu s
s = audioread("xpospi95.wav");

% Úloha 2
% Priemerovanie hodnôt z každých 16 vzoriek
i = 1;
j = 1;
while(i <= N)
  if(mean(s(i:i+16-1)) > 0)
    symbols(j) = 1;
  else
    symbols(j) = 0;
  endif
  i += 16; % Next symbol
  j += 1;
endwhile

% Kreslenie grafu
plot(linspace(0,0.02,320), s(1:320), "b");
title ("Signál a dekódované symboly");
xlabel("čas [ms]");
ylabel("symboly, s[n]");
hold on 
stem(linspace(0,0.02,20), symbols(1:20), "g");
hold off