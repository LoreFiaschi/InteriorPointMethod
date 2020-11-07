%R1 = [sqrt(2) 0 -1; 0 sqrt(3) 0; 1 0 sqrt(2)];
%R2 = [1 -1 0; 1 1 0; 0 0 sqrt(2)];
%R = 1/sqrt(6)*R2*R1;

%x = [0; 0; 3];

%y = R*x;

plot3([1 1], [1 2], [1 0]);
hold on
%grid on
plot3([1 2], [1 1], [1 0]);

%plot3([2 y(1)+2], [2 y(2)+2], [0 y(3)]);
%plot3([0 y(1)], [0 y(2)], [0 y(3)]);

R1 = [ sqrt(2/3) 0 -1/sqrt(3);
          0      1     0    ;
      1/sqrt(3) 0 sqrt(2/3)];
  
R2 = [ 1/sqrt(2) -1/sqrt(2) 0;
      1/sqrt(2) 1/sqrt(2) 0;
          0          0     1];
      
R = R2*R1;

%R = [ 1/sqrt(2) 1/sqrt(3) -1/sqrt(6);
%     -1/sqrt(2) 1/sqrt(3) -1/sqrt(6);
%      0         1/sqrt(3) -sqrt(2/3)];

L = zeros(3, 20);
z = linspace(0, 1, 20);
L(3,:) = z;

RL = R*L;

%plot3(L(1,:), L(2,:), L(3,:));
plot3(RL(1,:)+2, RL(2,:)+2, RL(3,:));

x = linspace(-0.5, 0.5, 20);
C = zeros(3, 40);
C(1,1:20) = x;
C(1,21:40) = x;
C(2,1:20) = sqrt(0.25-x.^2);
C(2,21:40) = -sqrt(0.25-x.^2);

RC = R*C;

RC(1:2,:) = RC(1:2,:)+2;

plot3(RC(1,:), RC(2,:), RC(3,:));

for i=1:3
    C(3,:) = C(3,:)+1/4*ones(1,40);
    RC = R*C;
    RC(1:2,:) = RC(1:2,:)+2;
    plot3(RC(1,:), RC(2,:), RC(3,:));
end

%gC = [2/3;1;-1/3]'*(C.^2)-2*sqrt(2/3)*C(1,:).*C(2,:)

grid on
xlabel('x')
ylabel('y')
zlabel('z')