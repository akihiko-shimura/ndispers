
from ndispers.media.crystals._KDP import KDP
import numpy as np

kdp = KDP()

print("Testing dndT...")
try:
    print(kdp.dndT(0.5, 0.1, 25, pol='o'))
except Exception as e:
    print(f"dndT failed: {e}")
