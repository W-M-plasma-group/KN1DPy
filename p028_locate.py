import numpy as np

def locate(table, value):
    """
    Finds the index of a value (or values) in a sorted table using np.searchsorted.
    Parameters:
        table (list or np.ndarray): A sorted list or array of numbers (ascending or descending).
        value (float, int, list, np.ndarray): Value(s) to search for in the table.        
    Returns:
        np.ndarray: An array of indices (integers) corresponding to the positions
                    where the values meet the conditions.
    """
    # Convert inputs to NumPy arrays if they are scalars or lists
    table = np.asarray(table)
    value = np.atleast_1d(value)  # Ensure `value` is an array
    
    # Determine if the table is in ascending or descending order
    asc = table[0] <= table[-1]
    
    if not asc:
        # If the table is in descending order, temporarily reverse it
        table = table[::-1]
    
    # Use np.searchsorted to find the indices
    indices = np.searchsorted(table, value, side='right' if asc else 'left') - 1
    # print(f"indices: {indices}")
    
    # Adjust indices for descending tables
    if not asc:
        indices = len(table) - indices - 1
    
    # Handle special cases: out-of-range values
    indices[value < table[0]] = -1  # Values less than the first element
    indices[value >= table[-1]] = len(table) - 1  # Values greater than or equal to the last element
    
    return indices

if __name__ == "__main__":
    table = [1, 3, 4, 6, 7, 8, 9]
    tbl_desc = [9, 8, 5, 4, 2, 1]
    
    print(locate(table, 2))       # Output: [0]
    print(locate(table, 5))       # Output: [2]
    print(locate(table, 3.4))     # Output: [1]
    print(locate(tbl_desc, 3.4))  # Output: [4]
    print(locate(table, [2, 5]))  # Output: [0, 2]
    print(locate(table, 0))       # Output: [-1]
    print(locate(table, 10))      # Output: [6]