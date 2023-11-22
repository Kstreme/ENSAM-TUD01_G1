using DifferentialEquations
using Plots

# Define the differential equation for y''''(x) = x
function ode!(du, u, p, x)
    du[1] = u[2]
    du[2] = u[3]
    du[3] = u[4]
    du[4] = p[1](x)  # Extract the function f(x) from the parameters
end

# Define the function f(x) = x
f(x) = sin(x)

# Set up the ODE problem for y''''(x) = x
a = 0.0  # Lower bound
b = 1.0  # Upper bound
u₀ = [1.0, 0.0, 0.0, 0.0]  # Initial values for y, y', y'', y'''
ode_prob = ODEProblem(ode!, u₀, (a, b), [f])

# Solve the ODE problem using Tsit5 solver
sol = solve(ode_prob, Tsit5())

# Extract the solution values for y(x)
x_values = range(a, b, length=100)
y_values = hcat(sol(x_values)...)  # Extract values for y, y', y'', y'''

# Plot the solution
plot(x_values, y_values[1, :], label="Solution", xlabel="x", ylabel="y(x)", legend=:top)
