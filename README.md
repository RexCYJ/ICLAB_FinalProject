# ICLAB_FinalProject
This is a part of the work for 2023 NTHUEE ICLab final project. Our group built a complete hybrid-crypto system for safe and reliable communication, which merged the symmetric (AES) and asymmetric (ECDH) encryption within a system.

We implemented and verified the algorithm from software to hardware, and, at the end, finished the ASIC design by 45-nm standard-cell library.

This repository collects my contribution to this work. I was responsible for implementing the encrytion and decryption function of <strong>ECDH (Elliptic Curve Diffie–Hellman key exchange)</strong>. 

Most of the algorithm and structure designs come from the reference paper below:

* [D. M. Schinianakis, A. P. Fournaris, H. E. Michail, A. P. Kakarountas and T. Stouraitis, "An RNS Implementation of an  $F_{p}$ Elliptic Curve Point Multiplier," in IEEE Transactions on Circuits and Systems I: Regular Papers, vol. 56, no. 6, pp. 1202-1213, June 2009, doi: 10.1109/TCSI.2008.2008507.][1]
* [J. Bucek and R. Lorencz, "Comparing Subtraction-Free and Traditional AMI," 2006 IEEE Design and Diagnostics of Electronic Circuits and systems, Prague, Czech Republic, 2006, pp. 95-97, doi: 10.1109/DDECS.2006.1649585.][2]

For more details of ECDH, here is my [notes][3]

[1]:<https://ieeexplore.ieee.org/document/4663678>
[2]:<https://ieeexplore.ieee.org/document/1649585>
[3]:<https://yujia-chen.notion.site/ECDH-77daa1fcd6b94cc298c2ff7273fe38af?pvs=4>
