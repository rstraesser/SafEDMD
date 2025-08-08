<!-- PROJECT SHIELDS -->
[![arXiv][arxiv-shield]][arxiv-url]
[![Webpage][webpage-shield-RS]][webpage-url-RS]

# SafEDMD: A Koopman-based data-driven controller design framework for nonlinear dynamical systems
`SafEDMD` stands for Stability- and certificate-oriented EDMD: a novel EDMD-based learning architecture which comes along with rigorous certificates, resulting in a reliable surrogate model generated in a data-driven fashion. Together with proportional error bounds, which vanish at the origin and are tailored to control tasks, it allows a certified controller design based on semi-definite programming.
## Installation
Download the `Matlab` files and install [[Yalmip]](https://yalmip.github.io/) as well as [[Mosek]](https://www.mosek.com/).

The provided code inspects two experiments:
* Nonlinear inverted pendulum: run `main_invertedPendulum`
* Nonlinear benchmark system: run `main_nonlinearInvariant`

## Reference
This repository contains an implementation of the ideas presented in the paper:

Strässer, R., Schaller, M., Worthmann, K., Berberich, J., & Allgöwer, F. "SafEDMD: A Koopman-based data-driven controller design framework for nonlinear dynamical systems", 2024, [[arxiv]](https://arxiv.org/abs/2402.03145)

---


If this software helped you with your research, please cite us.
```
@article{SafEDMD2024,
  title = {{SafEDMD}: A Koopman-based data-driven controller design framework for nonlinear dynamical systems},
  author = {Str{\"a}sser, Robin and Schaller, Manuel and Worthmann, Karl and Berberich, Julian and Allg{\"o}wer, Frank},
  year = {2024},
  journal={arXiv:2402.03145},
}

## Contact
🧑‍💻 Robin Strässer - [robin.straesser@ist.uni-stuttgart.de](mailto:robin.straesser@ist.uni-stuttgart.de)


[webpage-shield-RS]: https://img.shields.io/badge/Webpage-Robin%20Strässer-T?style=flat&logo=codementor&color=green
[webpage-url-RS]: https://www.ist.uni-stuttgart.de/institute/team/Straesser/
[arxiv-shield]: https://img.shields.io/badge/arXiv-2402.03145-t?style=flat&logo=arxiv&logoColor=white&color=red
[arxiv-url]: https://arxiv.org/abs/2402.03145


