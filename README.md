<a name="readme-top"></a>

[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
<!-- [![GPLv3 License][license-shield]][license-url] -->
<!-- [![LinkedIn][linkedin-shield]][linkedin-url] -->



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/zezhong-zhang/coredensity">
    <img src="logo.jpeg" alt="Logo" width="80" height="80">
  </a>

  <h3 align="center">Core Density</h3>

  <p align="center">
    Providing electron density for all element
    <br />
    <!-- <a href="https://coredensity.readthedocs.io/en/latest/"><strong>Explore the docs »</strong></a> -->
    <br />
    <br />
    <a href="https://github.com/zezhong-zhang/coredensity/tree/main/density.py">View Demo</a>
    ·
    <a href="https://github.com/zezhong-zhang/coredensity/issues">Report Bug</a>
    ·
    <a href="https://github.com/zezhong-zhang/coredensity/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

<!-- [![Product Name Screen Shot][product-screenshot]](https://github.com/zezhong-zhang/coredensity) -->
The electron density of an atom represents the probability distribution of orbital electrons as a function of distance from the nucleus. The accurate electron density is critical for scattering physics & materials characterization: the Fourier transform of the density is the X-ray scattering factor, and one can further deduce the electron scattering factor by using the Mott formula.

As the core electron of the heavy element can be very fast -- approaching the speed of light, one has to consider the relativistic effects. This project provides the electron density of each orbital for all the elements calculated by solving the Dirac equation.

The project aims to provide accurate scattering factors for electrons and X-rays, with the freedom to do any combination of orbitals. 

This repo is currently under development. 

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- GETTING STARTED -->
## Getting Started

To get code running locally, let's first create a conda environment.

```bash
conda create -n dense
pip install h5py numpy matplotlib
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ROADMAP -->
## Roadmap
### TO-DO list
- [x] Provide the density database solved from Dirac equation
- [x] Generate the X-ray scattering factor
- [x] Generate the electron scattering factor
- [ ] Combine the core electrons from Dirac and the rest valence electrons from DFT to consider the effect of bonding
- [ ] Generate the projected potential 
- [ ] Provide an interface to multislice simulation (aiming to be general for any packages)

See the [open issues](https://github.com/zezhong-zhang/coredensity/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions are what makes the open-source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the GPLv3 License. See `LICENSE` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Zezhong Zhang - zezhong.zhang@uantwerpen.be

Project Link: [https://github.com/zezhong-zhang/coredensity](https://github.com/zezhong-zhang/coredensity)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments
* To write....
* Readme page based on [Best-README-Template](https://github.com/othneildrew/Best-README-Template)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/zezhong-zhang/coredensity.svg?style=for-the-badge
[contributors-url]: https://github.com/zezhong-zhang/coredensity/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/zezhong-zhang/coredensity.svg?style=for-the-badge
[forks-url]: https://github.com/zezhong-zhang/coredensity/network/members
[stars-shield]: https://img.shields.io/github/stars/zezhong-zhang/coredensity.svg?style=for-the-badge
[stars-url]: https://github.com/zezhong-zhang/coredensity/stargazers
[issues-shield]: https://img.shields.io/github/issues/zezhong-zhang/coredensity.svg?style=for-the-badge
[issues-url]: https://github.com/zezhong-zhang/coredensity/issues
[license-shield]: https://img.shields.io/github/license/zezhong-zhang/coredensity.svg?style=for-the-badge
[license-url]: https://github.com/zezhong-zhang/coredensity/blob/main/LICENSE
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/zezhong-zhang-062a0838
[product-screenshot]: images/screenshot.png