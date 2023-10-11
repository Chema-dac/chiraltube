# chiraltube

Made by José María de Albornoz Caratozzolo & Felipe Cervantes Sodi.

If you use this program please cite it as such:   J. M. de Albornoz-Caratozzolo and F. Cervantes-Sodi, Nanoscale Adv., 2023, DOI: 10.1039/D3NA00301A

Chiraltube is a small python code that prints atomic coordinates for different nanotube and nanoribbon structures in special .xyz files for their usage in other visualization, simulation or calculation software.

Usage: `python3 ~/chiraltube.py <input file> [<output file>] [<options>]`

Input files should contain the unit cell of a specific 2D material in special .xyz format or in Quantum Espresso-style .in input files. Several examples of usable unit cells can be found in `~/chiraltube/UnitCellExamples`.

Use option `-h` / `-help` for more help.
