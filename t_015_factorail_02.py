from t_015_factorail_01 import shared_factorial

def calculate(n):
    # Incrementar el contador en cada llamada
    shared_factorial.call_count += 1
    if n == 0:
        return 1
    return n * calculate(n - 1)  # Llamada recursiva