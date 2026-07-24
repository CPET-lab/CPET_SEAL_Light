# CPET-SEAL Light 1.0

- **Noise Flooding**
- **Sparse Secret**


# 1. Implementation

### 1) Noise Flooding
- Implemented the core functions for bootstrapping.
- Source files:
  - [`noiseflooding.h`](/native/src/seal/noiseflooding.h)

### 2) Keygenerator & RLWE (Modified SEAL code)
- Implemented sparse ternary polynomial generation with a specified Hamming weight.
- Source files:
  - [`keygenerator.h`](/native/src/seal/keygenerator.h)
  - [`keygenerator.cpp`](/native/src/seal/keygenerator.cpp)
  - [`rlwe.h`](/native/src/seal/util/rlwe.h)
  - [`rlwe.cpp`](/native/src/seal/util/rlwe.cpp)
 
