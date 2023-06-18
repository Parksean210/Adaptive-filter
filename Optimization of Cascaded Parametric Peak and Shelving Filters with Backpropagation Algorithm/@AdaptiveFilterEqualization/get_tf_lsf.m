function tf_lsf = get_tf_lsf(obj, G, fc)
%METHOD1 이 메서드의 요약 설명 위치
%   자세한 설명 위치
V0 = 10^(G/20);
H0 = V0 - 1;
if G > 0 % boost
    a = (tan(pi*fc/obj.fs) - 1)/(tan(pi*fc/obj.fs) + 1);
else % cut
    a = (tan(pi*fc/obj.fs) - V0)/(tan(pi*fc/obj.fs) + V0);
end
A1 = obj.z1*[a; 1]./obj.z1*[1; a];
tf_lsf = 1 + H0/2*(1 + A1);
end
