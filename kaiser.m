%% creating a filter
% index no 170389M
A=3;B=8;C=9;

ap=0.03+(0.01*A); %maximum passband ripple
aa=45+B; % minimum stopband attenuation

omega_p1=(C*100)+400; % lower passband edge
omega_p2=(C*100)+950;  % upper passband edge
omega_a1=(C*100)+500; %lower stopband edge
omega_a2=(C*100)+800;  %upper stopband edge
omega_s=2*((C*100)+1300);  % sampling frequency

%finding cutoff frequencies
bt= min((omega_a1-omega_p1),(omega_p2-omega_a2));
omega_c1=omega_p1+(bt/2);
omega_c2=omega_p2-(bt/2);       

%finding delta
delta_p=((10^(0.05*ap))-1)/((10^(0.05*ap)) + 1);
delta_a = 10^(-0.05*aa);
delta=min(delta_p,delta_a);

%finding actual stopband attenuation
attenuation=-20*log10(delta);

%find alpha
if (attenuation<=21)
    alpha=0;
elseif (attenuation>21) && (attenuation<=50)
    alpha = 0.5842*(attenuation-21).^0.4 + 0.07886*(attenuation-21); 
else
    alpha = 0.1102*(attenuation-8.7);
end
%alpha=0;

%finding D
if (attenuation<=21)
    D=0.9222;
else
    D=(attenuation-7.95)/14.36;
end

%finding N
if mod(ceil(omega_s*D/bt)+1,2)==1      %ceil function
    N=ceil(omega_s*D/bt)+1;
elseif mod(ceil(omega_s*D/bt)+1,2)==0
    N=ceil(omega_s*D/bt)+2;
end

%find beta
n_t = (N-1)/2;  
n = -n_t:1:n_t; 
beta = alpha*((1-(2*n/(N-1)).^2).^0.5); 

%% kaiser window plotting
i_beta = 1; 
i_alpha = 1; 

limit=150;
for k = 1 : limit
    i_beta = i_beta + ((1/factorial(k))*(beta/2).^k) .^2; 
    i_alpha = i_alpha + ((1/factorial(k))*(alpha/2).^k).^2; 
end


window = i_beta./i_alpha ;
figure;
stem(n,window);
xlabel('n');
ylabel('Window[n]') ; 
title ( 'Windowing Function'); 

%% ideal impulse response filter
T=(2*pi)/omega_s;

n_negative = -n_t:1:-1;
h_negative = (1./(pi*n_negative)).*(sin(omega_c1*T.*n_negative)-sin(omega_c2*T.*n_negative)) ;

n_positive = 1:1:n_t;
h_positive = (1./(pi*n_positive)).*(sin(omega_c1*T.*n_positive)-sin(omega_c2*T.*n_positive));

n0=0;
h0 = 1+2*(omega_c1-omega_c2)/omega_s;

h = [h_negative,h0,h_positive];
n = [n_negative,n0,n_positive];

figure; 
stem(n,h);
xlabel('n');
title('ideal impulse response(Bandstop)');
ylabel('h[n]');

%% non-causal
h_filtered=h.*window;
figure;
stem(n,h_filtered);
xlabel('n');
ylabel('h[n]');
title('impulse response of a non causal filter');

%% causal
n_causal=[0:1:(n_t*2)];
figure;
stem(n_causal,h_filtered);
title('impulse response of a causal filter');
xlabel('n');
ylabel('h[n]');

%% frequency response
fvtool(h_filtered);
freqz(h_filtered);%obtaining the frequency response and corresponding frequencies

%%

%define frequencies of input signals%
omega_1=omega_c1/2;
omega_2=omega_c1+(omega_c2-omega_c1)/2;
omega_3=omega_c2+(omega_s/2-omega_c2)/2;
ns=0:1:300;
excitation= cos(omega_1*T.*ns)+cos(omega_2*T.*ns)+cos(omega_3*T.*ns);

%plot input signal- time domain
figure;
stem(ns,excitation);
title('Input signal - time domain');
xlabel('n');
ylabel('amplitude');
%%
%frequency domain - input signal
n_y_1=2^nextpow2(numel(ns));
y=fft(excitation,n_y_1)
func_x=omega_s/2*linspace(0,1,n_y_1/2+1);
figure;
plot(func_x,2*abs(y(1:n_y_1/2+1)));
title('input signal - frequency domain');
xlabel('frequency');
ylabel('X(f)');

%frequency domain - bandstop filter
func_h=fft(h_filtered,n_y_1) 
figure;
plot(func_x,2*abs(func_h(1:n_y_1/2+1)));
title('filter - frequency domain');
xlabel('frequency');
ylabel('X(f)');

%%output through filter
out=y.*func_h;
figure;
plot(func_x,2*abs(out(1:n_y_1/2+1)));
title('output - frequency domain');
xlabel('frequency');
ylabel('X(f)');


out_t=ifft(out);
figure;
stem(out_t(1:300));
title('output signal - time domain');
xlabel('n');
ylabel('amplitude)');
%%

%to check plot ideal response
%don't send the signal which is in bandstop frequency
%time domain
x_out=cos(omega_1*T.*ns)+cos(omega_3*T.*ns);
figure;
stem(ns,x_out);
title('Ideal filtered signal(desired) - time domain');
xlabel('n');
ylabel('amplitude');

%frequency domain
y=fft(x_out,n_y_1) 
func=omega_s/2*linspace(0,1,n_y_1/2+1);
figure;
plot(func,2*abs(y(1:n_y_1/2+1)));
title('ideal filtered signal(desired) - frequency domain');
xlabel('frequency');
ylabel('X(f)');

