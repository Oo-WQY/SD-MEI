# SD-MEI Bloch Equation Solver

## Overview
**SD-MEI (Spectral-Diagonalization-based Matrix Exponential Integration)** is a highly efficient algorithm for solving the Bloch equation, primarily used in forward modeling of surface nuclear magnetic resonance (sNMR). Compared to traditional numerical methods such as RK4, SD-MEI improves computational efficiency by one to two orders of magnitude, significantly accelerating the modeling process, especially for long-duration sequences like SSFP.

---

## File Descriptions

### 1. `F_norm_error_of_H.mlx`
- **Function**: Error analysis of the spectral decomposition of the system matrix **H**.
- **Description**: A key module for evaluating the error level of the SD-MEI algorithm. The **Frobenius norm** is used to quantify decomposition errors.

### 2. `SDMEI_Bloch_Mt.mlx`
- **Function**: Time-domain solution of the magnetization vector (*M*<sub>*t*</sub>).
- **Supported Sequences**: FID, SE, SSFP.
- **Usage**:  
  Modify the following parameter to select the sequence type:  
  `jlmrsData.measure.pulse_type = 'SSFP'`  
  Replace `'SSFP'` with `'FID'` or `'SE'` as needed.

### 3. `SDMEI_Bloch_B1_M.mlx`
- **Function**: Analysis of the magnetization vector's dependence on the excitation **B<sub>1</sub>** field.
- **Supported Sequences**: FID, SSFP.
- **Usage**:  
  Modify the following parameter to switch between sequence types:  
  `jlmrsData.measure.pulse_type = 'SSFP'`  
  Replace `'SSFP'` with `'FID'` as needed.

### 4. `SSFP_KernelLossFig.mlx`
- **Function**: Visualization of forward modeling kernel loss for SSFP sequences.

### 5. `SSFP_e0LossFig.mlx`
- **Function**: Visualization of the SSFP sequence signal (*e*<sub>0</sub>) loss.

---

## Usage Notes
- For sequence-related scripts (e.g., `SDMEI_Bloch_Mt.mlx` and `SDMEI_Bloch_B1_M.mlx`), you can quickly switch between simulation sequences by adjusting the **`pulse_type`** parameter.
- All files are in MATLAB Live Script format (`.mlx`), supporting interactive execution and result visualization.
