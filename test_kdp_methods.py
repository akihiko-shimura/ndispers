
from ndispers.media.crystals._KDP import KDP
import numpy as np

kdp = KDP()

print("Testing n...")
try:
    # This fails currently
    print(kdp.n(0.5, 0.1, 25, pol='e'))
except Exception as e:
    print(f"n failed: {e}")

print("Testing woa_theta...")
try:
    # This might fail if symbols mismatch args
    print(kdp.woa_theta(0.5, 0.1, 25, pol='e'))
except Exception as e:
    print(f"woa_theta failed: {e}")
