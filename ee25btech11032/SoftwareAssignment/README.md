#Image Compression Usig Randomized SVD and Jacobi Method

This project implements **image compression in C** using the **Randomized Singular Value Decomposition (rSVD)** algorithm, followed by the **One Sided Jacobi** algorithm.

---

## Goal

The goal of this project is to reuce the dimensionality of image data by decomposing the image matrix into low-rank approximations.It uses **Randomized SVD** for faster computation on large matrices and **Jacobi eigenvalue decomposition** to refine the results.

---

## Functional Overview

### This project is capable of
1. **Compressing Greyscale Iamges** (Channel = 1) to different degrees depending on different target rank (k') values.
2. **Converting a Coloured Image** (Channels = 3)to grey scale image and compress it further.
3. **Compressing a Coloured Image**

It also provides the flexibility to:
- Use **rSVD followed by Jacobi**, or  
- Perform compression **entirely using the Jacobi method**.

---

### Project Reports Compares
** k = 5 **
** k = 20 **
** k = 50 **
** k = 100 **

---
