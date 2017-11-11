/*
LAMMPS Script Maker - Bilipid Membrane with Hole, DPD Fluid (clark_bowman@brown.edu)
For LAMMPS version Aug. 10, 2015
*/

#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265359

double randDouble()
{
    return double(rand()) / double(RAND_MAX);
}

int main()
{
    srand(5000);
    using namespace std;

    double box_width, membrane_center, intra_spacing, box_length, temperature, v_lipid_bond, v_lipid_angle, v_lattice_density, v_timestep, hole_radius, finite_border;
    int target_lipids, tail_length, num_tails, num_heads, v_run_length, v_dump_interval, dump_tail_length;
    bool fill_fluid;

    // Read parameters from text file
    {
        ifstream file_in ("params_membrane_finite.txt");
        file_in >> box_width; // Y, Z widths of simulation box (LJ units)
        file_in >> box_length; // X length of simulation box
        file_in >> membrane_center; // Starting X position of membrane
        file_in >> intra_spacing; // Distance in between atoms of lipid tails
        file_in >> target_lipids; // Target number of lipids in the membrane
        file_in >> tail_length; // Length of each lipid tail
        file_in >> num_tails; // Number of hydrophobic tails per lipid
        file_in >> num_heads; // Number of hydrophilic heads per lipid
        file_in >> temperature; // LJ temperature of system
        file_in >> v_lipid_bond; // Harmonic bond strength in lipid
        file_in >> v_lipid_angle; // Harmonic angle strength in lipid tails
        file_in >> v_lattice_density; // Numerical density of fluid initialization
        file_in >> v_timestep; // LJ time per simulation step
        file_in >> v_run_length; // Number of steps to rain data-gathering simulation
        file_in >> v_dump_interval; // Dump visualization data every # steps
        file_in >> hole_radius; // Radius, in LJ units, of hole in center of membrane
        file_in >> fill_fluid; // Whether to fill hole with fluid (1) or not (0)
        file_in >> finite_border;
        file_in >> dump_tail_length;
    }

    double inter_spacing = box_width / sqrt(double(target_lipids / 2));
    int n_iter = int(ceil((box_width - 2. * finite_border) / inter_spacing));
    hole_radius = pow(hole_radius, 2);

    // Generate membrane initial condition
    {
        stringstream fileTemp;

        // Masses for each particle type
        fileTemp << "Masses" << endl << endl;
        fileTemp << "1 1.0" << endl << "2 1.0" << endl << "3 1.0" << endl << "4 1.0" << endl;

        // Loop to create atoms
        fileTemp << endl << "Atoms" << endl << endl;

        int atom_no = 1;
        int mol_no = 1;
        double ang, xpos, myx, myy, myz, yoff, zoff, ang_perp, ang_perp_inc;
        double ypos = finite_border;

        // Looping over all Y positions for lipids...
        for (int yctr = 0; yctr < n_iter && mol_no <= target_lipids; ypos += inter_spacing)
        {
            double zpos = finite_border;

            // Looping over all Z positions for lipids...
            for (int zctr = 0; zctr < n_iter && mol_no <= target_lipids; zpos += inter_spacing)
            {
                if (pow(ypos - box_width * 0.5, 2) + pow(zpos - box_width * 0.5, 2) > hole_radius)
                {
                    // Start with X position at the membrane center, then subtract off the length of a tail
                    xpos = membrane_center - (double(tail_length) + 0.5) * intra_spacing;
                    yoff = 0.;
                    zoff = 0.;

                    // For each head particle per lipid...
                    for (int hctr = 0; hctr < num_heads; hctr++)
                    {
                        // Create the head atom
                        fileTemp << atom_no + hctr << " " << mol_no << " 1 " << xpos << " " << ypos + yoff << " " << zpos + zoff << endl;

                        // If enough tails haven't been made yet...
                        if (hctr < num_tails)
                        {
                            // Initialize the tail location
                            myx = xpos;
                            myy = ypos + yoff;
                            myz = zpos + zoff;

                            // Create each tail atom, moving along the length of the tail in the X direction
                            for (int j = 0; j < tail_length; j++)
                            {
                                myx += intra_spacing;
                                fileTemp << atom_no + num_heads + hctr * tail_length + j << " " << mol_no << " " << 2 + 2 * int(j < dump_tail_length) << " " << myx << " " << myy << " " << myz << endl;
                            }
                        }

                        // Choose a random angle and offset the next head (if present) in that direction
                        ang = randDouble() * 2 * PI;
                        yoff += cos(ang) * intra_spacing;
                        zoff += sin(ang) * intra_spacing;
                    }

                    // Update atom count for number of atoms in one lipid, increment molecule count by 1
                    atom_no += num_heads + num_tails * tail_length;
                    mol_no++;

                    // Now do the same thing for the opposite layer of the membrane:
                    // Start with X position at the membrane center, then add the length of a tail
                    xpos = membrane_center + (double(tail_length) + 0.5) * intra_spacing;
                    yoff = 0.;
                    zoff = 0.;
                    for (int hctr = 0; hctr < num_heads; hctr++)
                    {
                        fileTemp << atom_no + hctr << " " << mol_no << " 1 " << xpos << " " << ypos + yoff << " " << zpos + zoff << endl;
                        if (hctr < num_tails)
                        {
                            myx = xpos;
                            myy = ypos + yoff;
                            myz = zpos + zoff;
                            for (int j = 0; j < tail_length; j++)
                            {
                                myx -= intra_spacing;
                                fileTemp << atom_no + num_heads + hctr * tail_length + j << " " << mol_no << " " << 2 + 2 * int(j < dump_tail_length) << " " << myx << " " << myy << " " << myz << endl;
                            }
                        }
                        ang = randDouble() * 2 * PI;
                        yoff += cos(ang) * intra_spacing;
                        zoff += sin(ang) * intra_spacing;
                    }
                    atom_no += num_heads + num_tails * tail_length;
                    mol_no++;
                }
                zctr++;
            }
            yctr++;
        }
        mol_no--;

        // Loop to create bonds
        fileTemp << endl << "Bonds" << endl << endl;
        int bond_no = 1;
        atom_no = 1;

        // For each lipid...
        for (int ctr = 0; ctr < mol_no; ctr++)
        {
            // The heads are bonded
            for (int i = 0; i < num_heads - 1; i++)
            {
                fileTemp << bond_no << " 1 " << atom_no + i << " " << atom_no + i + 1 << endl;
                bond_no++;
            }
            for (int i = 0; i < num_tails; i++)
            {
                // Each tail is bonded to a head
                fileTemp << bond_no << " 1 " << atom_no + i << " " << atom_no + num_heads + i * tail_length << endl;
                bond_no++;

                // Atoms within a tail are bonded
                for (int j = 0; j < tail_length - 1; j++)
                {
                    fileTemp << bond_no << " 1 " << atom_no + num_heads + i * tail_length + j << " " << atom_no + num_heads + i * tail_length + j + 1 << endl;
                    bond_no++;
                }
            }
            atom_no += (num_heads + num_tails * tail_length);
        }

        // Loop to create angles
        fileTemp << endl << "Angles" << endl << endl;
        int ang_no = 1;
        atom_no = 1;

        // For each lipid...
        for (int ctr = 0; ctr < mol_no; ctr++)
        {
            // Add angles along the length of each tail (no angle potential is imposed among heads, even if there are 3 or more)
            for (int i = 0; i < num_tails; i++)
            {
                fileTemp << ang_no << " 1 " << atom_no + i << " " << atom_no + num_heads + i * tail_length << " " << atom_no + num_heads + i * tail_length + 1 << endl;
                ang_no++;
                for (int j = 0; j < tail_length - 2; j++)
                {
                    fileTemp << ang_no << " 1 " << atom_no + num_heads + i * tail_length + j << " " << atom_no + num_heads + i * tail_length + j + 1 << " " << atom_no + num_heads + i * tail_length + j + 2 << endl;
                    ang_no++;
                }
            }
            atom_no += (num_heads + num_tails * tail_length);
        }
        ofstream fileOut ("membrane_finite.dat");
        fileOut.precision(5);
        fileOut << "# Membrane data file" << endl << endl;

        // Header with number of atoms, bonds, angles, dihedrals, and simulation information (dihedrals are not used here)
        fileOut << atom_no - 1 << " atoms" << endl;
        fileOut << bond_no - 1 << " bonds" << endl;
        fileOut << ang_no - 1 << " angles" << endl;
        fileOut << 0 << " dihedrals" << endl << endl;
        fileOut << 4 << " atom types" << endl << "1 bond types" << endl;
        fileOut << "1 angle types" << endl << "0 dihedral types" << endl << endl;
        fileOut << 0. << " " << box_length << " xlo xhi" << endl;
        fileOut << 0. << " " << box_width << " ylo yhi" << endl;
        fileOut << 0. << " " << box_width << " zlo zhi" << endl << endl;
        fileOut << fileTemp.rdbuf();
    }

    // Create LAMMPS script
    {
        ofstream fileOut ("membrane_finite.in");
        fileOut.precision(5);

        fileOut << "# LAMMPS Script - Bilipid Membrane with Hole, DPD Fluid (clark_bowman@brown.edu)" << endl;
        fileOut << "# LAMMPS version Aug. 10, 2015" << endl << endl << endl << endl << endl << endl;

        fileOut << "# SEC: This section initializes the simulation." << endl << endl;
        fileOut << "# 	Define units, atom style, log path, and neighbor settings; read configuration data for the membrane." << endl;
        fileOut << "# 	Final line communicates ghost data and is necessary for DPD parallelizing. On some versions of LAMMPS, use instead `communicate single vel yes'" << endl << endl;
        fileOut << "units lj" << endl;
        fileOut << "atom_style molecular" << endl;
        fileOut << "log membrane.log" << endl;
        fileOut << "read_data membrane_finite.dat" << endl;
        fileOut << "neighbor 0.3 bin" << endl;
        fileOut << "neigh_modify delay 3" << endl;
        fileOut << "comm_modify vel yes" << endl << endl;

        fileOut << "# 	Define bond, angle, pairwise interactions." << endl;
        fileOut << "# 	For pairwise interactions: type 1 is lipid head, 2 is lipid tail, 3 is water." << endl;
        fileOut << "#   43872 is a temperature seed; change for different randomness." << endl << endl;
        fileOut << "bond_style harmonic" << endl;
        fileOut << "bond_coeff 1 " << v_lipid_bond << " " << intra_spacing << endl;
        fileOut << "angle_style cosine/delta" << endl;
        fileOut << "angle_coeff 1 " << v_lipid_angle << " 180.0" << endl;
        fileOut << "dihedral_style none" << endl << endl;
        fileOut << "pair_style dpd " << temperature << " 1.0 43872" << endl;
        fileOut << "pair_coeff * * 25.0 4.5 1.0" << endl;
        fileOut << "pair_coeff 1 3 35.0 4.5 1.0" << endl;
        fileOut << "pair_coeff 2 3 75.0 20.0 1.0" << endl;
        fileOut << "pair_coeff 3 4 75.0 20.0 1.0" << endl;
        fileOut << "pair_coeff 1 2 50.0 9.0 1.0" << endl;
        fileOut << "pair_coeff 1 4 50.0 9.0 1.0" << endl;
        fileOut << "pair_coeff 2 2 15.0 4.5 1.0" << endl;
        fileOut << "pair_coeff 4 4 15.0 4.5 1.0" << endl;
        fileOut << "pair_coeff 2 4 15.0 4.5 1.0" << endl << endl;

        fileOut << "# SEC: This section initializes the geometry." << endl << endl;
        fileOut << "#   Define safe regions where fluid may be placed. This is the union of two regions outside the membrane." << endl;
        fileOut << "# 	The simulation box is " << box_length << " x " << box_width << " x " << box_width << "." << endl << endl;
        if (!fill_fluid)
        {
            fileOut << "region safe block " << membrane_center - double(tail_length + 1) * intra_spacing << " " << membrane_center + double(tail_length + 1) * intra_spacing << " " << finite_border << " " << box_width - finite_border << " " << finite_border << " " << box_width - finite_border << " units box side out" << endl;

        }
        else
        {
            fileOut << "region safe1 block " << membrane_center - double(tail_length + 1) * intra_spacing << " " << membrane_center + double(tail_length + 1) * intra_spacing << " " << finite_border << " " << box_width - finite_border << " " << finite_border << " " << box_width - finite_border << " units box side out" << endl;
            fileOut << "region safe2 cylinder x " << box_width * 0.5 << " " << box_width * 0.5 << " " << max(sqrt(hole_radius) * 0.95, sqrt(hole_radius) - 0.5) << " 1 " << box_length - 1 << " units box" << endl;
            fileOut << "region safe union 2 safe1 safe2" << endl;
        }
        fileOut << "lattice fcc " << v_lattice_density << endl;
        fileOut << "create_atoms 3 region safe" << endl << endl;

        fileOut << "# SEC: This section defines LAMMPS groups and computes that will be used in the simulation." << endl << endl;
        fileOut << "# 	Groups are named indicatively of their membership." << endl;
        fileOut << "# 	Computes include the x, y, and z positions of the center of mass of the membrane." << endl << endl;
        fileOut << "group fluid type 3" << endl;
        fileOut << "group membrane type 1 4" << endl << endl;

        fileOut << "# SEC: This section initializes the particle velocities." << endl << endl;
        fileOut << "velocity all create " << temperature << " 21456" << endl << endl;

        fileOut << "# SEC: This section defines fixes to impose forces in the simulation." << endl << endl;
        fileOut << "# 	NVE integration, with limit for initialization." << endl << endl;
        fileOut << "fix 1 all nve/limit 0.05" << endl << endl;

        fileOut << "# SEC: This section runs the simulation." << endl << endl;
        fileOut << "# 	Simulation timestep (LJ time units)." << endl;
        fileOut << "timestep " << v_timestep << endl << endl;
        fileOut << "# 	How often to output thermo data and first run phase." << endl;
        fileOut << "thermo 1000" << endl;
        fileOut << "run 100" << endl << endl;

        fileOut << "# 	Dump 1 is a fluid-less simulation dump at the specified interval." << endl << endl;
        fileOut << "dump 1 membrane atom " << v_dump_interval << " membrane_finite.lammpstrj" << endl;
        fileOut << "dump_modify 1 first yes" << endl << endl;

        fileOut << "# 	Run data-collecting simulation." << endl;
        fileOut << "run " << v_run_length << endl;
    }
    return 0;
}
