# SVD Image Compression 🖼️📉

This repository contains a university group project that implements an image compression tool using Singular Value Decomposition (SVD), developed entirely from scratch in Python. The system processes an input image, computes the SVD via custom mathematical algorithms, and reconstructs the image at various approximation ranks to demonstrate the relationship between data compression and visual fidelity.

## 🚀 Key Features
* **Zero-Dependency SVD:** The core decomposition algorithm relies exclusively on native Python structures and the `math` module, avoiding pre-built linear algebra libraries like NumPy.
* **Jacobi Rotations:** Implements an iterative Jacobi eigenvalue algorithm to compute the $U$, $V$, and $\Sigma$ matrices.
* **Custom Matrix Engine:** Features a full suite of manual matrix operations, including multiplication, transposition, scalar operations, and Householder transformations.
* **Image Processing Pipeline:** Includes custom functions to convert 3D RGB matrices to 2D grayscale matrices by averaging color channels, followed by downsampling for computational efficiency.
* **Rank-Based Reconstruction:** Generates approximations of the original image using truncated singular value ranks (e.g., $r=5$, $r=20$, $r=32$).

## 🛠️ Prerequisites
* **Python 3.x**
* `matplotlib` (used solely for reading the input image and plotting the final results)

Install the required visualization dependency via terminal:
```bash
pip install matplotlib
```

## ⚙️ Usage
1. Clone the repository and navigate to the project directory.
2. Ensure you have an input image named `picture.jpg` in the root folder.
3. Run the compression script:
```bash
python script_plot_image.py
```

## 📊 Outputs
Upon execution, the script generates a series of `matplotlib` figures:
1. **Grayscale Conversion:** The original image converted to a grayscale matrix.
2. **Rank Approximations:** Reconstructed images truncated at various ranks to show progressive visual improvements.
3. **Singular Values Plot:** A logarithmic plot detailing the decay of the singular values.
4. **Cumulative Variance:** A cumulative sum curve illustrating the variance captured by the retained singular values.

---
**Acknowledgments:** 
This project was coordinated by Lect. dr. Mihai Bucătaru at the Facultatea de Matematică și Informatică, Universitatea din București.