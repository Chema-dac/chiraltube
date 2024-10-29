# chiraltube

Made by José María de Albornoz Caratozzolo & Felipe Cervantes Sodi.

If you use this program please cite it as such:   J. M. de Albornoz-Caratozzolo and F. Cervantes-Sodi, Nanoscale Adv., 2023, DOI: 10.1039/D3NA00301A

Chiraltube is a small python code that prints atomic coordinates for different nanotube and nanoribbon structures in special .xyz files for their usage in other visualization, simulation or calculation software.

Usage: `python3 ~/chiraltube.py <input file> [<output file>] [<options>]`

Input files should contain the unit cell of a specific 2D material in special .xyz format or in Quantum Espresso-style .in input files. Several examples of usable unit cells can be found in `~/chiraltube/UnitCellExamples`.

Use option `-h` / `-help` for more help.

### Now supports Multi-Walled Nanotube generation!
If you use the option `-mw` you will be prompted for the number of layers to build and then the chiral indices.

Afterwards you will be prompted to select a layer-scaling option, these are explained here in more detail:

* Option 1: Doesn't transform any layer, returns the full MW NT as it was first created. (Recommended if all layers have the same chirality, and therefore have the same height, e.g. both are zigzag NTs)
* Option 2: Leaves biggest layer unchanged, all the other layers are repeated along the z-axis to surpass the biggest layer. Then they are trimmed so that they are all the same height. (Recommended if you don't care about periodicity of the full NT)
* Option 3: Leaves biggest layer unchanged, all other layers are scaled up along the z-axis to match the height of the biggest one. (Recommended if all NT heights are very close to each other)
* Option 4: Leaves biggest layer unchanged. Combination of option 2 and 3. First repeats the smaller layers along the z-axis to get as close as possible to the biggest height, then scales up or down to match the remaining difference. (Recommended, although it might scale too much in certain cases, losing physical plausibility)
* Option 5: Only works with two layers for now. Looks for a combination of numbers that, when repeating each layer by its corresponding number, makes both layers get close enough for scaling to be practical. (Recommended in most cases, preserves periodicity while staying close to the original stucture. Might not find a supercell, which makes it default to option 4)
