# A neutron source term for OpenMC simulations of the UoB High-Flux Accelerator-Driven Neutron Facility (HF-ADNeF)
**Author:** Max Conroy



> [!IMPORTANT]
>This source term has been designed to be used with the HF-ADNeF target-room geometry available on the [UoB Nuclear Group github](https://github.com/bhamnuclear/hfadnef-openmc-geometry). As such, neutrons are initialised at $z=-4.51$ cm, and on a 5 cm radius disc in $x,y$. The code can be edited to produce different positional distributions, for example a point source.


### Overview
This respository contains the code for a [Parameterised Compiled Source](https://docs.openmc.org/en/latest/usersguide/settings.html) for use in OpenMC simulations of HF-ADNeF at the University of Birmingham. The source generates neutrons produced via the $^7Li(p,n)$ reaction via a Monte Carlo method. This code can be used as a starting source term for neutrons produced via the <sup>7</sup>Li(p,n)<sup>7</sup>Be reaction, which is the production method at HF-ADneF.

For full details of how the code works, please see *paper in progress*.

In brief:
- Neutrons are produced uniformly in x and y over a 5-cm-radius disk, centred at x = 0, y = 0.
- The z position of the neutron production depends on energy loss calculated from SRIM and the total cross section of the interaction. Near threshold, this uses the description given by Lee and Zhou [^1] and above this energy, experimental data from EXFOR are used.
- The emission vector of the neutrons is sampled from the differential cross section data given by Liskien and Paulsen [^2].
- The neutrons are given a weight such that tally results are per mC of proton current.

### Installation
OpenMC will need to be installed on your system. This can be done via the instructions given here: https://docs.openmc.org/en/latest/quickinstall.html. I have found it most reliabl to manually build OpenMC from source using CMake.

With OpenMC installed the CompiledSource can now be compiled agains it. The /src folder contains the .cpp file which samples neutron starting conditions from physical principles and experimentally measured data, which are stored in /src/data.

To use the code, first compile it via:

```cmake -S src -B build```

Then:

```cmake --build build```

For more information on using compiled sources, see the OpenMC documentation: https://docs.openmc.org/en/stable/usersguide/settings.html.

### Usage:
You can then use the libsource.so file as a compiled starting source in an OpenMC simulation. 

A ```CompiledSource``` can be initialised with the following line:
```
settings.source = openmc.CompiledSource('<path_to_source>/build/libsource.so', Ep)
```
The parameter `Ep` is the energy of the proton beam, which should be passed as a string. The maximum energy that can be used is Ep = 2.6 MeV.

If you wish to add some spread to the energy of the incident proton beam, you will need to edit the `n_source.cpp` file. Simply uncomment/comment the relevant sections within the main sampling loop and re-compile.

**See the Jupyter notebook in the examples folder for a simple use case.**

> [!TIP]
> If you use this code in your work, please reference it accordingly. I will be writing this up as part of my PhD and will update this github when that happens. Please contact me for more information at m.j.conroy@pgr.bham.ac.uk.

### References:
[^1]: C. L. Lee and X. L. Zhou. “Thick target neutron yields for the 7Li(p,n)7Be reaction near threshold”. In: Nuclear Instruments and Methods in Physics Research Section B: Beam Interactions with Materials and Atoms 152.1 (Apr. 1999), pp. 1–11. issn: 0168-583X. doi: 10.1016/S0168-583X(99)00026-9.
[^2]: Horst Liskien and Arno Paulsen. “Neutron production cross sections and energies for the reactions 7Li(p,n)7Be and 7Li(p,n)7Be*”. In: Atomic Data and Nuclear Data Tables 15.1 (Jan. 1975), pp. 57–84. issn: 0092-640X. doi: 10.1016/0092-640X(75)90004-2.
