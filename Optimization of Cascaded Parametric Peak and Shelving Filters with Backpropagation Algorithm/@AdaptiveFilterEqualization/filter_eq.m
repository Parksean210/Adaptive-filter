function filter_eq(obj)
%FILTER_EQ 이 함수의 요약 설명 위치
%   자세한 설명 위치
if isempty(obj.sos)
    error('SOS must be calculated to use this function')
end
num_filter = length(obj.type);
num_sample = length(obj.input);
state = zeros(num_filter, 2);
obj.target = zeros(num_sample, 1);
for i = 1:num_sample
    x = obj.input(i);
    for j = 1:num_filter
        [x, state(j, :)] = filter(obj.sos(j,1:3), obj.sos(j,4:6), x, state(j,:));
    end
    obj.target(i) = x;
end
end

