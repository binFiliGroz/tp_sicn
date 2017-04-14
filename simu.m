function teb = simu (rapport_Eb_N0, N)
T = 0.000001;
F = 16;
L = 8;
alpha = 0.5;

bn = round(rand(1, N));

ak=2*bn.-1;

t_a = 0 : T : (N-1)*T;

s_t = zeros(1, N*F);
s_t (1:F:length(s_t)) = ak * F;
t_s = 0 : T/F : (N*F-1)*T/F;

t_filtre = [0 : T/F : L*T - (T/F)];

h_t = gen_filters3('srrc', t_filtre, T, F, L, alpha);

x_t = conv(s_t, h_t, "same");
L_xt = length(x_t);
t_x = 0 : T/F : L_xt*T/F-T/F;

Px_t = mean(x_t.^2);

rapport_Eb_N0_lin = 10^(rapport_Eb_N0/10);
sigma_n = sqrt(Px_t*F/2*1/rapport_Eb_N0_lin);
n_t = sigma_n*randn(1, length(x_t));

r_t = n_t + x_t;

h_r = fliplr(h_t);
p_t = conv(h_t, fliplr(h_t));

t2=T/F :T/F :2*(L*T)-T/F;

y_t = conv(r_t, h_r, "same");
taille_visu = 2;
begin_offset = length(h_t);
end_offset = 2*length(h_t);

t0 = 0; 
y_k = downsample(y_t, F, t0);

b_n = sign(y_k);

erreurs = sum(xor((ak+1)/2, (b_n+1)/2));
teb = erreurs/N;
endfunction
