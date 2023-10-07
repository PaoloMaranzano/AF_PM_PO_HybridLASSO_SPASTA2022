function [Harm,Harm_names] = FourierBasisMat(x,K,period)

%%%%% Frequencies
freq = 2*pi*[1:K]/period

%%%%% Harmonics
Harm = zeros(length(x),2*K);
Harm_names = strings(1,3);
for i = 1:length(freq)
    Harm(:,(i*2-1)) = sin(freq(i)*x);
    Harm_names((i*2-1)) = ['S' num2str(i) '_' num2str(period)];
    Harm(:,(i*2)) = cos(freq(i)*x);
    Harm_names((i*2)) = ['C' num2str(i) '_' num2str(period)];
end

end


