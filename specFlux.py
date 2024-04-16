#%%
from sympy import symbols, Eq, solve

def calculate_missing_value(s2=None, s1=None, freq2=None, freq1=None, alpha=None):
    """
    Calculate the missing value in the equation (s2/s1) = (freq2/freq1)**alpha.

    Parameters:
        s2 (float): The value of s2 (default is None).
        s1 (float): The value of s1 (default is None).
        freq2 (float): The value of freq2 (default is None).
        freq1 (float): The value of freq1 (default is None).
        alpha (float): The value of alpha (default is None).

    Returns:
        The calculated missing value or a message indicating no solution.
    """
    
    # Define symbols for s2, s1, freq2, freq1, and alpha
    s2_sym, s1_sym, freq2_sym, freq1_sym, alpha_sym = symbols('s2 s1 freq2 freq1 alpha')
    
    # Create the equation (s2/s1) = (freq2/freq1)**alpha
    equation = Eq(s2_sym / s1_sym, (freq2_sym / freq1_sym) ** alpha_sym)
    
    # Prepare a dictionary to hold the known values
    known_values = {}
    
    # Check for the provided values and add them to the dictionary
    if s2 is not None:
        known_values[s2_sym] = s2
    if s1 is not None:
        known_values[s1_sym] = s1
    if freq2 is not None:
        known_values[freq2_sym] = freq2
    if freq1 is not None:
        known_values[freq1_sym] = freq1
    if alpha is not None:
        known_values[alpha_sym] = alpha
    
    # Solve the equation for the unknown value
    solution = solve(equation.subs(known_values), dict=True)
    
    # Check if a solution was found
    if not solution:
        return "No solution found for the given inputs."
    
    # Extract the missing value from the solution
    missing_variable = None
    
    for key in [s2_sym, s1_sym, freq2_sym, freq1_sym, alpha_sym]:
        if key not in known_values:
            missing_variable = key
            break
    
    # If a missing variable is found, return its value from the solution
    if missing_variable is not None:
        return solution[0].get(missing_variable, "No solution for the missing variable.")
    
    return "All values provided; no missing variable to solve for."

# Example usage:
s2 = 0.2324  # Given value
s1 = None  # Missing value
freq2 = 1.051  # Given value
freq1 = 0.835  # Given value
alpha = -0.4  # Given value

# Calculate the missing value
result = calculate_missing_value(s2=s2, s1=s1, freq2=freq2, freq1=freq1, alpha=alpha)
print(f"The result is {result}")

# %%
