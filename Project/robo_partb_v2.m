%% *** Robot (kinematic) model parameters *** 
    clear all; 
    close all; 
    %%Link values in cm
    l0 = 2.0;  l1 = 4.0;  
    l2 = 5.0;  l3 =3.0; 


    %% *** sampling period *** 
    %% *** for the robot motion, kinematic simulation: 
    dt = 0.005; %i.e. 1 msec)   

    %% *** Create (or load from file) reference signals *** 
    %% *** DESIRED MOTION PROFILE - TASK SPACE *** 
    Tf=10.0; 	% 10sec duration of motion A->B
    t=linspace(0,2*Tf,2*Tf/dt);  
    
    %A :start point , B:final point  --> desired task-space trajectory  
    A = [3.0;-2.0;8.0] ; B = [-4.0 ; 3.0; 8.0];  
    %Trajectory : linear segment (x0,y0)-->(x1,y1); Time duration: Tf; 
    disp('Initialising Desired Task-Space Trajectory (Motion Profile) ...'); %% 
    disp(' ');
    %Calculation of the polynomials and trajectory construction
    [x_efAB,y_efAB,ux_efAB,uy_efAB] = trajectory(A,B,Tf,dt);
    [x_efBA,y_efBA,ux_efBA,uy_efBA] = trajectory(B,A,Tf,dt);
    z_ef = A(3).*ones(length(t),1); %position of end effector on z-axis
    x_ef = cat(1,x_efAB',x_efBA'); % position of end effector on x-axis 
    y_ef = cat(1,y_efAB',y_efBA'); % position of end effector on y-axis
    ux_ef = cat(1,ux_efAB',ux_efBA'); % linear velocity of end effector on the x-axis 
    uy_ef = cat(1,uy_efAB',uy_efBA'); % linear velocity of end effector on the y-axis
    
   
%% Ploting The Position and the velocity of the End effector through time     
fig1 = figure;

subplot(2,2,1); 
plot(t,x_ef); 
ylabel('X Position of End Effector (cm)'); 
xlabel('Time (sec)');

subplot(2,2,2); 
plot(t,y_ef); 
ylabel('Y Position of End Effector (cm)'); 
xlabel('Time (sec)');

subplot(2,2,3); 
plot(t,ux_ef); 
ylabel('ux (cm/sec)'); 
xlabel('Time (sec)');

subplot(2,2,4); 
plot(t,uy_ef); 
ylabel('uy (cm/sec)'); 
xlabel('Time (sec)');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %% ****** KINEMATIC SIMULATION - Main loop ****** 
    disp('Kinematic Simulation ...'); %% 
    disp(' '); %%  
    
    %% ***** INVESRE KINEMATICS  -->  DESIRED MOTION - JOINT SPACE ***** 
    tmp =(x_ef.^2-l1.^2+(z_ef-l0).^2);
    p = real(sqrt(tmp));
    nom = (x_ef.^2)+(y_ef.^2)+((z_ef-l0).^2)-(l1.^2)-(l2.^2)-(l3.^2);
    denom = (2*l2*l3);
    q3 = acos(nom./denom);
    q1=atan2(z_ef-l0,x_ef)-atan2(l1,+(p)) ;
    q2=atan2(y_ef,+(p))-atan2(real(sin(q3(:))*l3),real(l2+cos(q3(:))*l3)); 
    
    
    qd1 = zeros(length(t),1); qd2=qd1; qd3=qd1;
for k = 1:length(t)
    [j,T] = Jinv(l0,l1,l2,l3,q1(k),q2(k),q3(k));
    qd1(k) = j(1,1)*ux_ef(k) + j(1,2)*uy_ef(k);
    qd2(k) = j(2,1)*ux_ef(k) + j(2,2)*uy_ef(k);
    qd3(k) = j(3,1)*ux_ef(k) + j(3,2)*uy_ef(k);
end

fig2 = figure;

subplot(2,3,1); 
plot(t,q1);
ylabel('q1 (rad)'); 
xlabel('time (sec)');

subplot(2,3,2); 
plot(t,q2); 
ylabel('q2 (rad)'); 
xlabel('time (sec)');

subplot(2,3,3); 
plot(t,q3);
ylabel('q3 (rad)'); 
xlabel('time (sec)');

subplot(2,3,4); 
plot(t,qd1);
ylabel('qd1 (rad)'); 
xlabel('time (sec)');

subplot(2,3,5); 
plot(t,qd2); 
ylabel('qd2 (rad)'); 
xlabel('time (sec)');

subplot(2,3,6); 
plot(t,qd3);
ylabel('qd3 (rad)'); 
xlabel('time (sec)');

    %% ***** FORWARD KINEMATICS  JOINT MOTION -->  CARTESIAN POSITIONS ***** 
c1 = cos(q1); c2 = cos(q2); c3 = cos(q3); c23 = cos(q2+q3);
s1 = sin(q1); s2 = sin(q2); s3 = sin(q3); s23 = sin(q2+q3);
    
x0=zeros(size(q1)); y0=zeros(size(q1)); z0=zeros(size(q1)); % Cartesian Position of base "joint"
x1=zeros(size(q1)); y1=zeros(size(q1)); z1=ones(size(q1)).*l0; % Cartesian Position of first joint 
x2=-l1.*sin(q1); y2=zeros(size(q1)); z2= l0+l1*cos(q1); % Cartesian Position of seconcd joint
x3=c1.*c2.*l2-l1.*s1;  y3=s2.*l2; z3=l0+c1.*l1+l2.*s1.*c2; % Cartesian Position of thrid joint 
xe=c1.*c23.*l3+c1.*c2.*l2-s1.*l1; % x-coordinate of end effector
ye=s23.*l3+s2.*l2; % y-coordiante of end effector 
ze=s1.*l3.*c23+s1.*c2.*l2+c1.*l1+l0; % z-coordinate of end effector 



fig3= figure; 
axis([-5 5 -5 5 0 10])                                              
axis on 
hold on 
grid on
title('Stick simulation of trajectory')
xlabel('x (cm)'); 
ylabel('y (cm)'); 
zlabel ('z (cm)');
dtk=100;                                                              



for tk=1:dtk:length(t) 
   pause(0.2);                                                           
   plot3([x0],[y0],[z0],'s');  
   plot3([x0(tk),x1(tk)],[y0(tk),y1(tk)], [z0(tk),z1(tk)]);
   plot3([x1],[y1],[z1],'o');  
   plot3([x1(tk),x2(tk)],[y1(tk),y2(tk)], [z1(tk),z2(tk)]);
   plot3([x2(tk)],[y2(tk)], [z2(tk)],'o');  
   plot3([x2(tk),x3(tk)],[y2(tk),y3(tk)],[z2(tk),z3(tk)]);	
   plot3([x3(tk)],[y3(tk)],[z3(tk)],'o');  
   plot3([x3(tk),xe(tk)],[y3(tk),ye(tk)],[z3(tk),ze(tk)]);
   plot3([xe(tk)],[ye(tk)],[ze(tk)],'o'); 
end
plot3(xe,ye,ze,'rs'); 

 
