function [A,B,C,D,z] = theta2abcd(theta, g)
    Xu = theta(1); Xq = theta(2);
    Mu = theta(3); Mq = theta(4);
    Xd = theta(5); Md = theta(6);
    A = [Xu, Xq, -g;
         Mu, Mq,  0;
         0,   1,  0];
    B = [Xd; Md; 0];
    C = [0 1 0];
    D = 0;
    if nargout == 5
        sys = ss(A,B,C,D);
        z = tzero(sys);
    end
end
