using DifferentialEquations
using Plots

function f!(du, u, p, x)
    du[1] = u[2]
    du[2] = u[3]
    du[3] = u[4]
    du[4] = x  # f(x) = x

    # No additional differential equation needed

end

u0 = [0.0, 0.0, 0.0, 0.0]  # Initial conditions
xspan = (0.0, 1.0)  # Define the x span for the solution

# Enforce u(x=0) = 0, u'(x=0) = 0, u''(x=1) = 0, and u'''(x=1) = 0 as boundary conditions
bc = (u, p, x) -> begin
    [u[1],     # u(x=0) - 0 = 0
     u[2],     # u'(x=0) - 0 = 0
     u[3] - 0, # u''(x=1) - 0 = 0
     u[4] - 0] # u'''(x=1) - 0 = 0
end

prob = ODEProblem(f!, u0, xspan, bc=bc)

# Solve the ODE problem
sol = solve(prob)

# Extract the solution at specific points
x_values = range(0, stop=1, length=100)  # Adjust as needed
u_values = [sol(x)[1] for x in x_values]
u_prime_values = [sol(x)[2] for x in x_values]
u_double_prime_values = [sol(x)[3] for x in x_values]
u_triple_prime_values = [sol(x)[4] for x in x_values]

# Create plots
plot(x_values, u_values, label="u")
plot!(x_values, u_prime_values, label="u'")
plot!(x_values, u_double_prime_values, label="u''")
plot!(x_values, u_triple_prime_values, label="u'''")
xlabel!("x")
ylabel!("Function values")
title!("Function and Derivatives")

# Display the plot
display(Plot)
