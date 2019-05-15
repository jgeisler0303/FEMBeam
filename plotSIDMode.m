function plotSIDMode(sid, mode)

for i= 1:length(sid.frame)
    X(:, i)= sid.frame(i).origin.M0;
    U(:, i)= sid.frame(i).Phi.M0(:, mode);
    R(:, i)= sid.frame(i).Psi.M0(:, mode);  
end

[~, extentIdx]= max(sum(X, 2));

subplot(6, 1, 1)
plot(X(extentIdx, :), U(1, :))
ylabel('x deflection')
grid on

subplot(6, 1, 2)
plot(X(extentIdx, :), U(2, :))
ylabel('y deflection')
grid on

subplot(6, 1, 3)
plot(X(extentIdx, :), U(3, :))
ylabel('z deflection')
grid on

subplot(6, 1, 4)
plot(X(extentIdx, :), R(1, :))
ylabel('x rot deflection')
grid on

subplot(6, 1, 5)
plot(X(extentIdx, :), R(2, :))
ylabel('y rot deflection')
grid on

subplot(6, 1, 6)
plot(X(extentIdx, :), R(3, :))
ylabel('z rot deflection')
grid on
