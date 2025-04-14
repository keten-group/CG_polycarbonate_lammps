# CG_polycarbonate
Materials for using Keten Group's coarse-grained polycarbonate model in LAMMPS

Last updated 04/14/2025 by Heather White

Source publications
- Wenjie Xia et al., Energy renormalization for coarse-graining polymers having different segmental structures. Science Advances (2019). DOI:10.1126/sciadv.aav4683

Related publications
- Nitin Hansoge et al., Effect of polymer chemistry on chain conformations in hairy nanoparitcle assemblies. ACS Macro Letters (2019). DOI: 10.1021/acsmacrolett.9b00526
- Amirhadi Alesadi et al., Understanding the role of cohesive interaction in mechanical behavior of a glass polymer. ACS Macromolecules (2020). DOI: 10.1021/acs.macromol.0c00067
- Nitin Hansoge et al., Universal relation for effective interaction between polymer-grafted nanoparticles. ACS Macromolecules (2021). DOI: 10.1021/acs.macromol.0c02600
- Jie Yang et al., Understanding the mechanical and viscoelastic properties of graphene reinforced polycarbonate nanocomposites using coarse-grained molecular dynamics simulations. Computational Materials Science (2021). DOI: 10.1016/j.commatsci.2021.110339

Polycarbonate is coarse-grained in a 3-bead representation as shown in the figure below. In previously developed materials, it is not uncommon for the beads 'A', 'B', and 'C' to be referred to as 'P', 'C', and 'O', respectively, so take care in keeping track of what naming convention is being used. For the purposes of this documentation, the 'ABC' representation will be used exclusively to avoid confusion. 

![image](https://github.com/keten-group/CG_polycarbonate/assets/55956937/ac579c5d-4e85-4da8-b5ab-d4349ed5bf8b)

Masses
- A: 76.0 g/mol
- B: 42.0 g/mol
- C: 60.0 g/mol

Bonds
- Harmonic style: https://docs.lammps.org/bond_harmonic.html
- Parameters
  - A-B / B-A: $K$ = 40.4 kcal/&#197;<sup>2</sup> , $r_0$ = 2.897 &#197;
  - A-C / C-A: $K$ = 112.3 kcal/&#197;<sup>2</sup>, $r_0$ = 3.425 &#197;

Angles
- Tabulated style: https://docs.lammps.org/angle_table.html
- Files:
  - ABA_angle_avg_potential.table
  - BAC_angle_avg_potential.table
  - ACA_angle_avg_potential.table
- Note: the tabulated potentials are fairly irregular and are not well approximated by harmonic functions.

Dihedrals
- Tabulated style: https://docs.lammps.org/dihedral_table.html
- File: ABAC_dihedral_avg_potential.table
- IMPORTANT: see notes on dihedral potential below

**Notes on the dihedral potential**: The PC molecules exhibit such high linearity that occasionally a dihedral angle is undefined. This produces extremely high forces that causes crashes with a "bond atoms missing" or similar error. This is not reproducable on a restart. The contribution of the dihedral potentials to the overall system is negligible and we recommend that these terms are neglected. 



