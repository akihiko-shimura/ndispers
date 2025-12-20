
from ndispers.media.crystals._KDP import KDP
import numpy as np

kdp = KDP()
try:
    print("Calling kdp.n(0.5, 0.1, 'e', pol='e')")
    print(kdp.n(0.5, 0.1, 'e', pol='e'))
except Exception as e:
    print(f"Error 3: {e}")
