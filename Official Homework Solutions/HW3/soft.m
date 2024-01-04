function z = soft(y,T)

y(abs(y) <= T) = 0;
y(y > T) = y(y > T) - T;
y(y < -T) = y(y < -T) + T;

z = y;