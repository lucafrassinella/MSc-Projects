function [A, B, C, D] = modelClass_SO(Xu, Xq, Mu, Mq, Xd, Md, ~)
% Single Output Model Class (output: q)

g = 9.81;
A = [Xu Xq -g; Mu Mq 0; 0 1 0];
B = [Xd Md 0]';
C = [0 1 0];
D = 0;

end