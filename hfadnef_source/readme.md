## HF-ADNeF Neutron Starting Source
This code can be used as a starting source term for neutrons produced via the <sup>7</sup>Li(p,n)<sup>7</sup>Be reaction, which is the production method at HF-ADneF.

For full details of how the code works, please see *document to be uploaded*.

In brief:
- Neutrons are produced uniformly in x and y over a 5-cm-radius disk, centred at y = -42 cm. This is to align with the proton beamline in CAD models of the facility.
- The z position of the neutron production depends on energy loss calculated from SRIM and the total cross section of the interaction. Near threshold, this uses the description given by Lee and Zhou [^1] and above this energy, experimental data from EXFOR are used.
- The emission vector of the neutrons is sampled from the differential cross section data given by Liskien and Paulsen [^2].
- The neutrons are given a weight such that tally results are per mC of proton current.

### Installation:
The src folder contains the .cpp file which samples neutron starting conditions from physical principles and experimentally measured data.

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

If you wish to *add an energy width* to the incident proton beam, you will need to edit the `n_source.cpp` file. Simply uncomment/comment the relevant sections within the main sampling loop.

**See the Jupyter notebook in the examples folder for a simple use case.**

> [!Caution]
> This code is a work in progress. There may be bugs and/or inaccuracies. Please get in contact if you notice anything which is not behaving as expected.

### References:
[^1]: C. L. Lee and X. L. Zhou. “Thick target neutron yields for the 7Li(p,n)7Be reaction near threshold”. In: Nuclear Instruments and Methods in Physics Research Section B: Beam Interactions with Materials and Atoms 152.1 (Apr. 1999), pp. 1–11. issn: 0168-583X. doi: 10.1016/S0168-583X(99)00026-9.
[^2]: Horst Liskien and Arno Paulsen. “Neutron production cross sections and energies for the reactions 7Li(p,n)7Be and 7Li(p,n)7Be*”. In: Atomic Data and Nuclear Data Tables 15.1 (Jan. 1975), pp. 57–84. issn: 0092-640X. doi: 10.1016/0092-640X(75)90004-2.
