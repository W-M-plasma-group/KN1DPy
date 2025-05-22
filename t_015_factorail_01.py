class FactorialShared:
    def __init__(self):
        self.call_count = 0  # Contador de llamadas recursivas

# Instancia global única
shared_factorial = FactorialShared()
