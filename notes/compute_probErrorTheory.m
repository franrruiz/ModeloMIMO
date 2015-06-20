% 2,4,6,8,10 QAM
flagEbNocte = 1;

s_e = [4 4 4 4];
s_p = [8 24 56 120];
s_c = [4 36 196 900];
d = [0.6325 0.3086 0.1524 0.0766];

s_t = s_e+s_p+s_c;

SNR_dB = 14;
s2y = 10^(-SNR_dB/10);

s2uni = s2y/2;
s2uni = s2uni./([4 6 8 10]/4);
stduni = sqrt(s2uni/20);

pe_e = 1-qfunct(-d./(2*stduni)).^2;
pe_p = 1-qfunct(-d./(2*stduni)).*(qfunct(-d./(2*stduni))-qfunct(d./(2*stduni)));
pe_c = 1-(qfunct(-d./(2*stduni))-qfunct(d./(2*stduni))).^2;

p_err = (pe_e.*s_e+pe_p.*s_p+pe_c.*s_c)./s_t;

p_err

