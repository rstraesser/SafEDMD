# SafEDMD: A certified learning architecture tailored to data-driven control of nonlinear dynamical systems
`SafEDMD` stands for Stability- and certificate-oriented EDMD: a novel EDMD-based learning architecture which comes along with rigorous certificates, resulting in a reliable surrogate model generated in a data-driven fashion. Together with proportional error bounds, which vanish at the origin and are tailored to control tasks, it allows a certified controller design based on semi-definite programming.
## Installation
Download the `Matlab` files and install [[Yalmip]](https://yalmip.github.io/) as well as [[Mosek]](https://www.mosek.com/).
## References
This repository contains an implementation of the ideas presented in the paper:

Strässer, R., Schaller, M., Worthmann, K., Berberich, J., & Allgöwer, F. "SafEDMD: A certified learning architecture tailored to data-driven control of nonlinear dynamical systems", 2024, [[arxiv]](https://arxiv.org/abs/2402.03145)

---


If this software helped you with your research, please cite us.
```
@article{SafEDMD2024,
  title = {{SafEDMD}: A certified learning architecture tailored to data-driven control of nonlinear dynamical systems},
  author = {Str{\"a}sser, Robin and Schaller, Manuel and Worthmann, Karl and Berberich, Julian and Allg{\"o}wer, Frank},
  year = {2024},
  journal={arXiv:2402.03145},
}
