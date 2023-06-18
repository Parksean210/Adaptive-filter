function [tf_pf, num, den] = get_tf_pf(obj, G, fb, fc)
%METHOD1 이 메서드의 요약 설명 위치
%   자세한 설명 위치
V0 = 10^(G/20);
H0 = V0 - 1;
d = -cos(2*pi*fc/obj.fs);
if G > 0 % boost
    a = (tan(pi*fb/obj.fs) - 1)/(tan(pi*fb/obj.fs) + 1);
else % cut
    a = (tan(pi*fb/obj.fs) - V0)/(tan(pi*fb/obj.fs) + V0);
end
num = [1+(1+a)*H0/2; d*(1-a); -a-(1+a)*H0/2];
den = [1; d*(1-a); -a];
tf_pf = obj.z2*num./(obj.z2*den);
end
