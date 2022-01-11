function [res,T] = Jinv(l0,l1,l2,l3,q1,q2,q3)

    % FORWARD KINEMATIC ANALYSIS
    A0h = [1,0,0,0 ; 0,0,-1,0 ; 0,1,0,l0 ; 0,0,0,1];
    Ah1 = [cos(q1),0,-sin(q1),0 ; sin(q1),0,cos(q1),0 ; 0,-1,0,0 ; 0,0,0,1];
    A12 = [cos(q2),-sin(q2),0,l2*cos(q2) ; sin(q2),cos(q2),0,l2*sin(q2) ; 0,0,1,l1 ; 0,0,0,1];
    A2E = [cos(q3),-sin(q3),0,l3*cos(q3) ; sin(q3),cos(q3),0,l3*sin(q3) ; 0,0,1,0 ; 0,0,0,1];
    A01 = A0h*Ah1; A02 = A01*A12; T = A02*A2E;
    % DIFFERENTIAL ANALYSIS + DIFFERENTIAL INVERSE ANALYSIS 
    b0h = A0h(1:3,3); b1 = A01(1:3,3); b2 = A02(1:3,3);
    p0h = A0h(1:3,4); p1 = A01(1:3,4); p2 = A02(1:3,4); pe = T(1:3,4);
    J = [cross(b0h,pe-p0h),cross(b1,pe-p1),cross(b2,pe-p2) ; b0h,b1,b2]; res = inv(J(1:3,:));
end