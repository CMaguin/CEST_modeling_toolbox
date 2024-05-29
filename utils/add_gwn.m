function s_noisy=add_gwn(s, i)
%This function adds gaussian white noise (i dB) to the signal s.

    sigma=sqrt(i);
    
    
    s_noisy=s + sigma* randn(size(s));