clear all;
close all;

N = 2048;
T = 0.000001;
F = 16;
L = 8;
alpha = 0.5;
rapport_Eb_N0 = 2;

# Question 1

bn = round(rand(1, N));

Mbn=mean(bn);
Vbn=var(bn);

# Question 2

ak=2*bn.-1;

Mak=mean(ak);
Vak=var(ak);
Pak=mean(ak.^2);

# Question 3

t_a = 0 : T : (N-1)*T;

figure(1);
stem(t_a, ak);
axis([0 (N-1)*T -1.2 1.2]);
title('Realisation de la suite {ak}');
xlabel('Temps[s]');
grid;
print("graphs/q3.png");

# Question 4

figure(2);
plot(real(ak), imag(ak),'x');
axis([-2 2 -2 2]);
title('Constellation');
xlabel('Re(ak)');
ylabel('Im(ak)');
print("graphs/q4.png");

# Question 5

s_t = zeros(1, N*F);
s_t (1:F:length(s_t)) = ak * F;
t_s = 0 : T/F : (N*F-1)*T/F;

figure(3);
stem(t_s, s_t); 
axis([0 (N-1)*T -1.2*F 1.2*F]);
title('Expansion');
xlabel('Temps[s]');
print("graphs/q5.png");

# Question 6

t_filtre = [0 : T/F : L*T - (T/F)];

h_t = gen_filters3('nrz', t_filtre, T, F, L, alpha);

figure(4);
#plot(t_filtre, h_t);
#title('Reponse NRZ');
#xlabel('Temps[s]');
#print("graphs/q6_nrz.png");
spectrum2(h_t, T/F, 'log');
print("graphs/q6_nrz_bis.png");

h_t = gen_filters3('rz', t_filtre, T, F, L, alpha);

figure(5);
#plot(t_filtre, h_t);
#title('Reponse RZ');
#xlabel('Temps[s]');
#print("graphs/q6_rz.png");
spectrum2(h_t, T/F, 'log');
print("graphs/q6_rz_bis.png");


h_t = gen_filters3('srrc', t_filtre, T, F, L, alpha);

figure(6);
#plot(t_filtre, h_t);
#title('Reponse SRRC');
#xlabel('Temps[s]');
#print("graphs/q6_srrc.png");
spectrum2(h_t, T/F, 'log');
print("graphs/q6_srrc_bis.png");

# Question 7

x_t = conv(s_t, h_t, "same");
L_xt = length(x_t);
t_x = 0 : T/F : L_xt*T/F-T/F;

figure(7);
subplot(211);
plot(t_s, s_t);
grid;
axis([0 N*T -1.5 1.5]);
title('signal avant filtre de mise en forme');
xlabel('Temps[s]');
ylabel('Amplitude [V]');
subplot(212);
plot(t_x, x_t, 'r');
grid;
axis([0 N*T -1.5 1.5]);
title('signal apres filtre de mise en forme he');
xlabel('Temps [s]');
ylabel('Amplitude [V]');
print("graphs/q7.png");

# Question 8

Px_t = mean(x_t.^2);

# Question 9

#figure(8);
#DSP_xt = spectrum2(xcorr(x_t), T/F, 'log');
#print("graphs/q9.png");

# Quesion 11

rapport_Eb_N0_lin = 10^(rapport_Eb_N0/10);
sigma_n = sqrt(Px_t*F/2*1/rapport_Eb_N0_lin);
n_t = sigma_n*randn(1, length(x_t));

r_t = n_t + x_t;

# Question 13

h_r = fliplr(h_t);
p_t = conv(h_t, fliplr(h_t));

t2=T/F :T/F :2*(L*T)-T/F;
figure(9);
plot(t2, p_t);
grid;
print("graphs/q13.png");

#plot(roots(p_t), "x");

# Question 14

y_t = conv(r_t, h_r, "same");
taille_visu = 2;
begin_offset = length(h_t);
end_offset = 2*length(h_t);

figure(10);
eyepattern(y_t,T,F,taille_visu,1,1);
print("graphs/q14a.png");

figure(11);
eyepattern(y_t,T,F,taille_visu,begin_offset,end_offset);
print("graphs/q14b.png");

# Question 18

t0 = 0; 
y_k = downsample(y_t, F, t0);
figure(11);
plot(real(y_k), imag(y_k),'x');
print("graphs/q19.png");

# Question 19

b_n = sign(y_k);
figure(12);
hold on;
stem(t_a, b_n);
plot(t_s, y_t, t_s, x_t);
plot(t_a, ak, "+", t_a, b_n, "x");
hold off;
print("graphs/q18.png");

erreurs = sum(xor((ak+1)/2, (b_n+1)/2));
teb = erreurs/N;
