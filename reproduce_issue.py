
from ndispers.media.crystals._KDP import KDP
import numpy as np

kdp = KDP()
try:
    print("Calling kdp.n(0.5, 0.1, pol='e')")
    print(kdp.n(0.5, 0.1, pol='e'))
except Exception as e:
    print(f"Error 1: {e}")

try:
    print("Calling kdp.n(0.5, 0.1, 'e')")
    print(kdp.n(0.5, 0.1, 'e'))
except Exception as e:
    print(f"Error 2: {e}")
