function y = h_stub(x, u)
% Predict measurement y = j
V = u(1);
y = 100 * sqrt(V);  % matches your current Plant Newton toy behavior
end