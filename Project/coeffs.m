function Result=coeffs(start_point,end_point,t0 , tf)
    %Equations % 
    %î1(t)=a0+a1t+a2t^2+a3t^3+a4t^4
    %î2(t)=b0+b1t
    %î3(t)=c0+c1t+c2t^2+c3t^3+c4t^4
    %x1(t)=î1(t) ,x2(t)=î2(t),x3(t)=î3(t)'
    
    syms a0 a1 a2 a3 a4 b0 b1 c0 c1 c2 c3 c4 ;
    Delta=(tf-t0)*0.05;
    DDelta=2*Delta;
    td=t0+DDelta;
    
    eq1=a0+a1*t0+a2*t0^2+a3*t0^3+a4*t0^4-start_point==0;                                  %Position Continuity at t0
    eq2=a1+2*a2*t0 + 3*a3*t0^2 + 4*a4*t0^3  == 0;                                  %Velocity Continuity at t0
    eq3=a0 + a1*td + a2*td^2 + a3*td^3 + a4*td^4  - (b0 + b1*td) == 0;     %Position Continuity at td
    eq4 = a1 + 2*a2*td + 3*a3*td^2 + 4*a4*td^3  -b1 == 0;                        %Velocity Continuity at td
    eq5 = 2*a2 + 6*a3*td + 12*a4*td^2== 0;                                             %Acceleration Continuity at td
    eq6=c0+c1*tf+c2*tf^2+c3*tf^3+c4*tf^4-end_point==0;                                      %Position Continuity at tf
    eq7 = c1 + 2*c2*tf + 3*c3*tf^2 + 4*c4*tf^3  == 0;                                %Velocity Continuity at tf
    eq8 = 2*c2 + 6*c3*tf + 12*c4*tf^2  == 0;                                            %Acceleration Continuity at tf
    eq9 = c0 + c1*(tf-DDelta) + c2*(tf-DDelta)^2 + c3*(tf-DDelta)^3 + c4*(tf-DDelta)^4  - (b0 + b1*(tf-DDelta)) == 0;         %Position Continuity at tf-DDelta
    eq10 = c1 + 2*c2*(tf-DDelta) + 3*c3*(tf-DDelta)^2 + 4*c4*(tf-DDelta)^3  - b1 == 0;                                                    %Velocity Continuity at tf-DDelta
    eq11= b0+b1*td - (a0 + a1*t0 + a2*t0^2 + a3*t0^3 + a4*t0^4)  - (Delta*b1) == 0;                                   %Ö-1
    eq12 = (c0 + c1*tf + c2*tf^2 + c3*tf^3 + c4*tf^4)  - (b0+b1*(tf-DDelta)) - (Delta*b1) == 0;                       %Ö-3

    equations=[ eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10 eq11 eq12];
    variables=[a4 a3 a2 a1 a0 b1 b0 c4 c3 c2 c1 c0];

    [A,b] = equationsToMatrix(equations,variables);
    Result=vpa(linsolve(A,b));
end