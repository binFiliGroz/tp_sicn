clear all;
close all;

N = 20;
T = 0.000001;
F = 16;
L = 8;
alpha = 0.5;
rapport_Eb_N0 = 10;

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
stem(t_a, ak, "filled");
title('Realisation de la suite {ak}');
xlabel('Temps[s]');
grid;
print("graphs/q3.png");

# Question 4

plot(real(ak), imag(ak),'x');
axis([-2 2 -2 2]);
title('Constellation');
xlabel('Re(ak)');
ylabel('Im(ak)');
#print("graphs/q4.png");

# Question 5

s_t = zeros(1, N*F);
s_t (1:F:length(s_t)) = ak * 16;
t_s = 0 : T : (N*F-1)*T;
plot(t_s, s_t); 
title('Expansion');
xlabel('Temps[s]');
print("graphs/q5.png");

# Question 6

t_filtre = [0 : T/F : L*T - (T/F)];
h_t = gen_filters3('srrc', t_filtre, T, F, L, alpha);
plot(t_filtre, h_t);
title('Reponse');
xlabel('Temps[s]');
print("graphs/q6.png");


# Question 7

x_t = conv(s_t,h_t, "same");
L_xt = length(x_t);
t_x = 0 : T : (L_xt-1)*T;
#subplot(211);
#plot(t_s, s_t);
#grid;
#axis([0 L_xt*T/F -1.5 1.5]);
#title('signal avant filtre de mise en forme');
#xlabel('Temps[s]');
#ylabel('Amplitude [V]');
#subplot(212);
#plot(t_x, x_t, 'r');
#grid;
#axis([0 L_xt*T/F -1.5 1.5]);
#title('signal apres filtre de mise en forme he');
#xlabel('Temps [s]');
#ylabel('Amplitude [V]');

# Question 8

Px_t = mean(x_t.^2);

# Question 9

#DSP_xt = spectrum2(xcorr(x_t), T/F, 'log');

# Quesion 11

rapport_Eb_N0_lin = 10^(rapport_Eb_N0/10);
sigma_n = sqrt(Px_t*F/2*1/rapport_Eb_N0_lin);
n_t = sigma_n*randn(1, length(x_t));

r_t = n_t + x_t;

# Question 13

h_r = fliplr(h_t);
p_t=conv(h_t, fliplr(h_t));

t2=T/F :T/F :2*(L*T)-T/F;
#plot(t2, p_t);
#grid;

#plot(roots(p_t), "x");

# Question 14

y_t = conv(r_t, h_r, "same");
taille_visu = 2;
begin_offset = length(h_t);
end_offset = 2*length(h_t);
#eyepattern(y_t,T,F,taille_visu,1,1);
#eyepattern(y_t,T,F,taille_visu,begin_offset,end_offset);

# Question 18

t0 = 0; 
y_k = downsample(y_t, F, t0);
#plot(real(y_k), imag(y_k),'x');

# Question 19

b_n = sign(y_k);
#stem(t_a*16, b_n);
#hold on;
#plot(t_s, y_t, t_s, x_t);
#hold off;
#plot(t_a, ak, "+", t_a, b_n, "x");

erreurs = sum(xor((ak+1)/2, (b_n+1)/2));
teb = erreurs/N;
