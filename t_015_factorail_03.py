from t_015_factorail_01 import shared_factorial
from t_015_factorail_02 import calculate

def main():
    n = 5
    # Reiniciar el contador antes de cada cálculo
    shared_factorial.call_count = 0
    result = calculate(n)
    print(f"Factorial de {n}: {result}")
    print(f"Número de llamadas recursivas: {shared_factorial.call_count}")

if __name__ == "__main__":
    main()